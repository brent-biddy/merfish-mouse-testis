include { samplesFrom; publishedPath } from './samplesheet'

process CLUSTER_SPATIALDATA_GPU {
    tag "${sample}"

    label 'gpu'

    publishDir { "${params.outdir}/${sample}/cluster_spatialdata_gpu" }, mode: 'copy'

    // The store is carried, not staged: this step reads only the table.
    input:
    tuple val(sample), val(zarr_store), path(table)
    path 'timer.py'

    output:
    tuple val(sample), val(zarr_store), path("${sample}.cluster_spatialdata_gpu.h5ad"), emit: sdata
    path "${sample}.cluster_spatialdata_gpu.timing.tsv", emit: timings

    script:
    """
    cluster_spatialdata_gpu.py --sample ${sample} --table_path ${table} --outdir .
    """

    stub:
    """
    touch ${sample}.cluster_spatialdata_gpu.h5ad
    touch ${sample}.cluster_spatialdata_gpu.timing.tsv
    """
}

workflow cluster_spatialdata_gpu {
    take:
    // a samplesheet, or tuple(sample, zarr_store, table) per sample
    input

    main:
    def ch_tables = samplesFrom(input, ['sample', 'zarr_store', 'table_path'])

    CLUSTER_SPATIALDATA_GPU(ch_tables, file("${projectDir}/bin/timer.py"))

    CLUSTER_SPATIALDATA_GPU.out.sdata
        .map { sample, zarr_store, table ->
            "${sample},${publishedPath(input, sample, zarr_store)}," +
            "${params.outdir}/${sample}/cluster_spatialdata_gpu/${table.name}"
        }
        .collectFile(name: 'cluster_spatialdata_gpu_samplesheet.csv', storeDir: params.outdir,
                     seed: 'sample,zarr_store,table_path', newLine: true, sort: true)

    emit:
    sdata = CLUSTER_SPATIALDATA_GPU.out.sdata // tuple(sample, zarr_store, table)
}
