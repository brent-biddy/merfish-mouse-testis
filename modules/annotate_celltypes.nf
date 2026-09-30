include { samplesFrom; publishedPath } from './samplesheet'

process ANNOTATE_CELLTYPES {
    tag "${sample}"

    publishDir { "${params.outdir}/${sample}/annotate_celltypes" }, mode: 'copy'

    // The store is carried, not staged: this step reads only the table.
    input:
    tuple val(sample), val(zarr_store), path(table)
    path reference
    path 'timer.py'

    output:
    tuple val(sample), val(zarr_store), path("${sample}.annotate_celltypes.h5ad"), emit: sdata
    path "${sample}.gene_overlap.tsv", emit: gene_overlap
    path "${sample}.annotate_celltypes.timing.tsv", emit: timings

    script:
    """
    annotate_celltypes.py --sample ${sample} --table_path ${table} --reference ${reference} --outdir .
    """

    stub:
    """
    touch ${sample}.annotate_celltypes.h5ad
    touch ${sample}.gene_overlap.tsv
    touch ${sample}.annotate_celltypes.timing.tsv
    """
}

workflow annotate_celltypes {
    take:
    // a samplesheet, or tuple(sample, zarr_store, table) per sample
    input

    main:
    def reference = file(params.reference
        ?: "${projectDir}/assets/reference/shami_human_testis_centroids.csv.gz")

    def ch_tables = samplesFrom(input, ['sample', 'zarr_store', 'table_path'])

    ANNOTATE_CELLTYPES(ch_tables, reference, file("${projectDir}/bin/timer.py"))

    ANNOTATE_CELLTYPES.out.sdata
        .map { sample, zarr_store, table ->
            "${sample},${publishedPath(input, sample, zarr_store)}," +
            "${params.outdir}/${sample}/annotate_celltypes/${table.name}"
        }
        .collectFile(name: 'annotate_celltypes_samplesheet.csv', storeDir: params.outdir,
                     seed: 'sample,zarr_store,table_path', newLine: true, sort: true)

    emit:
    sdata = ANNOTATE_CELLTYPES.out.sdata // tuple(sample, zarr_store, table)
}
