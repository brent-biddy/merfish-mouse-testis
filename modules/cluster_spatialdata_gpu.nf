include { samplesFrom } from './samplesheet'

process CLUSTER_SPATIALDATA_GPU {
    tag "${sample}"

    label 'gpu'

    publishDir { "${params.outdir}/${sample}/cluster_spatialdata_gpu" }, mode: 'copy'

    input:
    tuple val(sample), path(zarr)
    path 'timer.py'

    output:
    tuple val(sample), path("${sample}.cluster_spatialdata_gpu.zarr"), emit: zarr
    path "${sample}.cluster_spatialdata_gpu.timing.tsv", emit: timings

    script:
    """
    cluster_spatialdata_gpu.py --sample ${sample} --path ${zarr} --outdir .
    """

    stub:
    """
    mkdir -p ${sample}.cluster_spatialdata_gpu.zarr
    touch ${sample}.cluster_spatialdata_gpu.timing.tsv
    """
}

workflow cluster_spatialdata_gpu {
    take:
    // a samplesheet, or tuple(sample, zarr) per sample
    input

    main:
    def ch_zarrs = samplesFrom(input, ['sample', 'path'])

    CLUSTER_SPATIALDATA_GPU(ch_zarrs, file("${projectDir}/bin/timer.py"))

    CLUSTER_SPATIALDATA_GPU.out.zarr
        .map { sample, zarr ->
            "${sample},${params.outdir}/${sample}/cluster_spatialdata_gpu/${zarr.name}"
        }
        .collectFile(name: 'cluster_spatialdata_gpu_samplesheet.csv', storeDir: params.outdir,
                     seed: 'sample,path', newLine: true, sort: true)

    emit:
    zarr = CLUSTER_SPATIALDATA_GPU.out.zarr
}
