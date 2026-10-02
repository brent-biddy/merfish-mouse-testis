include { samplesFrom } from './samplesheet'

process PREP_CELLPOSE_VPT {
    tag "${sample}"

    publishDir { "${params.outdir}/${sample}/prep_cellpose_vpt" }, mode: 'copy'

    input:
    tuple val(sample), path(region_dir), path(cellpose_dir), val(label_offset)
    path 'timer.py'

    output:
    tuple val(sample), path("*.{csv,parquet,json}"), emit: vpt_files
    path "${sample}.prep_cellpose_vpt.timing.tsv", emit: timings

    script:
    """
    prep_cellpose_vpt.py --sample ${sample} --path ${region_dir} \\
        --cellpose_path ${cellpose_dir} --label_offset=${label_offset} --outdir .
    """

    stub:
    """
    touch cellpose_cell_by_gene.csv cellpose_cell_metadata.csv cellpose_micron_space.parquet
    touch detected_transcripts.csv prep_cellpose_vpt.json
    touch ${sample}.prep_cellpose_vpt.timing.tsv
    """
}

workflow prep_cellpose_vpt {
    take:
    // a samplesheet, or tuple(sample, region dir, cellpose dir) per sample
    input

    main:
    def ch_segmentations = samplesFrom(input, ['sample', 'path', 'cellpose_path'])

    // label_offset is optional and not a path, so samplesFrom cannot take it: a blank or
    // missing column means merge.py's raster is the whole mosaic. The step checks it either way.
    def ch_offsets = input instanceof Path || input instanceof String
        ? channel.fromPath(input).splitCsv(header: true, quote: '"')
            .map { row -> tuple(row.sample, row.label_offset ?: '0,0') }
        : ch_segmentations.map { row -> tuple(row[0], '0,0') }

    PREP_CELLPOSE_VPT(ch_segmentations.join(ch_offsets), file("${projectDir}/bin/timer.py"))

    // Rejoined rather than carried through the process, which never reads it: the region
    // staged into the task is a work dir, not where the caller pointed.
    PREP_CELLPOSE_VPT.out.vpt_files
        .join(ch_segmentations.map { sample, region_dir, cellpose_dir -> tuple(sample, "${region_dir}") })
        .map { sample, files, region_path ->
            "${sample},${region_path},${params.outdir}/${sample}/prep_cellpose_vpt"
        }
        .collectFile(name: 'prep_cellpose_vpt_samplesheet.csv', storeDir: params.outdir,
                     seed: 'sample,path,vpt_path', newLine: true, sort: true)

    emit:
    vpt_files = PREP_CELLPOSE_VPT.out.vpt_files
}
