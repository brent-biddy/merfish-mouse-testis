include { samplesFrom } from './samplesheet'

process PREP_CELLPOSE_VPT {
    tag "${sample}"

    publishDir { "${params.outdir}/${sample}/prep_cellpose_vpt" }, mode: 'copy'

    input:
    tuple val(sample), path(region_dir), path(cellpose_dir)
    path 'timer.py'

    output:
    tuple val(sample), path("*.{csv,parquet}"), emit: vpt_files
    path "${sample}.prep_cellpose_vpt.timing.tsv", emit: timings

    script:
    """
    prep_cellpose_vpt.py --sample ${sample} --path ${region_dir} \\
        --cellpose_path ${cellpose_dir} --outdir .
    """

    stub:
    """
    touch cellpose_cell_by_gene.csv cellpose_cell_metadata.csv cellpose_micron_space.parquet
    touch detected_transcripts.csv
    touch ${sample}.prep_cellpose_vpt.timing.tsv
    """
}

workflow prep_cellpose_vpt {
    take:
    // a samplesheet, or tuple(sample, region dir, cellpose dir) per sample
    input

    main:
    def ch_segmentations = samplesFrom(input, ['sample', 'path', 'cellpose_path'])

    PREP_CELLPOSE_VPT(ch_segmentations, file("${projectDir}/bin/timer.py"))

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
