include { samplesFrom } from './samplesheet'

process ANNOTATE_CELLTYPES {
    tag "${sample}"

    publishDir { "${params.outdir}/${sample}/annotate_celltypes" }, mode: 'copy'

    input:
    tuple val(sample), path(zarr)
    path reference
    path 'timer.py'

    output:
    tuple val(sample), path("${sample}.annotate_celltypes.zarr"), emit: zarr
    path "${sample}.gene_overlap.tsv", emit: gene_overlap
    path "${sample}.annotate_celltypes.timing.tsv", emit: timings

    script:
    """
    annotate_celltypes.py --sample ${sample} --path ${zarr} --reference ${reference} --outdir .
    """

    stub:
    """
    mkdir -p ${sample}.annotate_celltypes.zarr
    touch ${sample}.gene_overlap.tsv
    touch ${sample}.annotate_celltypes.timing.tsv
    """
}

workflow annotate_celltypes {
    take:
    // a samplesheet, or tuple(sample, zarr) per sample
    input

    main:
    def reference = file(params.reference
        ?: "${projectDir}/assets/reference/shami_human_testis_centroids.csv.gz")

    def ch_zarrs = samplesFrom(input, ['sample', 'path'])

    ANNOTATE_CELLTYPES(ch_zarrs, reference, file("${projectDir}/bin/timer.py"))

    ANNOTATE_CELLTYPES.out.zarr
        .map { sample, zarr ->
            "${sample},${params.outdir}/${sample}/annotate_celltypes/${zarr.name}"
        }
        .collectFile(name: 'annotate_celltypes_samplesheet.csv', storeDir: params.outdir,
                     seed: 'sample,path', newLine: true, sort: true)

    emit:
    zarr = ANNOTATE_CELLTYPES.out.zarr
}
