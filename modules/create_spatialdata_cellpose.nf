include { samplesFrom } from './samplesheet'

process CREATE_SPATIALDATA_CELLPOSE {
    tag "${sample}"

    publishDir { "${params.outdir}/${sample}/create_spatialdata_cellpose" }, mode: 'copy'

    input:
    tuple val(sample), path(region_dir), path(vpt_dir)
    path 'timer.py'

    output:
    tuple val(sample), path("${sample}.create_spatialdata_cellpose.zarr"), emit: zarr
    path "${sample}.create_spatialdata_cellpose.timing.tsv", emit: timings

    script:
    """
    create_spatialdata_cellpose.py --sample ${sample} --path ${region_dir} \\
        --vpt_path ${vpt_dir} --outdir .
    """

    stub:
    """
    mkdir -p ${sample}.create_spatialdata_cellpose.zarr
    touch ${sample}.create_spatialdata_cellpose.timing.tsv
    """
}

workflow create_spatialdata_cellpose {
    take:
    // a samplesheet, or tuple(sample, region dir, vpt dir) per sample
    input

    main:
    def ch_region_dirs = samplesFrom(input, ['sample', 'path', 'vpt_path'])

    CREATE_SPATIALDATA_CELLPOSE(ch_region_dirs, file("${projectDir}/bin/timer.py"))

    CREATE_SPATIALDATA_CELLPOSE.out.zarr
        .map { sample, zarr ->
            "${sample},${params.outdir}/${sample}/create_spatialdata_cellpose/${zarr.name}"
        }
        .collectFile(name: 'create_spatialdata_cellpose_samplesheet.csv', storeDir: params.outdir,
                     seed: 'sample,path', newLine: true, sort: true)

    emit:
    zarr = CREATE_SPATIALDATA_CELLPOSE.out.zarr
}
