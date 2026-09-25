include { samplesFrom } from './samplesheet'

process CREATE_SPATIALDATA {
    tag "${sample}"

    publishDir { "${params.outdir}/${sample}/create_spatialdata" }, mode: 'copy'

    input:
    tuple val(sample), path(region_dir)
    path 'timer.py'

    output:
    tuple val(sample), path("${sample}.create_spatialdata.zarr"), emit: zarr
    path "${sample}.create_spatialdata.timing.tsv", emit: timings

    script:
    """
    create_spatialdata.py --sample ${sample} --path ${region_dir} --outdir .
    """

    stub:
    """
    mkdir -p ${sample}.create_spatialdata.zarr
    touch ${sample}.create_spatialdata.timing.tsv
    """
}

workflow create_spatialdata {
    take:
    // a samplesheet, or tuple(sample, region dir) per sample
    input

    main:
    def ch_region_dirs = samplesFrom(input, ['sample', 'path'])

    CREATE_SPATIALDATA(ch_region_dirs, file("${projectDir}/bin/timer.py"))

    CREATE_SPATIALDATA.out.zarr
        .map { sample, zarr -> "${sample},${params.outdir}/${sample}/create_spatialdata/${zarr.name}" }
        .collectFile(name: 'create_spatialdata_samplesheet.csv', storeDir: params.outdir,
                     seed: 'sample,path', newLine: true, sort: true)

    emit:
    zarr = CREATE_SPATIALDATA.out.zarr
}
