include { samplesFrom } from './samplesheet'

process CREATE_SPATIALDATA {
    tag "${sample}"

    publishDir { "${params.outdir}/${sample}/create_spatialdata" }, mode: 'copy'

    input:
    tuple val(sample), path(region_dir)
    path 'timer.py'

    output:
    tuple val(sample), path("${sample}.create_spatialdata.zarr"),
          path("${sample}.create_spatialdata.h5ad"), emit: sdata
    path "${sample}.create_spatialdata.timing.tsv", emit: timings

    script:
    """
    create_spatialdata.py --sample ${sample} --path ${region_dir} --outdir .
    """

    stub:
    """
    mkdir -p ${sample}.create_spatialdata.zarr
    touch ${sample}.create_spatialdata.h5ad
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

    CREATE_SPATIALDATA.out.sdata
        .map { sample, zarr_store, table ->
            def published = "${params.outdir}/${sample}/create_spatialdata"
            "${sample},${published}/${zarr_store.name},${published}/${table.name}"
        }
        .collectFile(name: 'create_spatialdata_samplesheet.csv', storeDir: params.outdir,
                     seed: 'sample,zarr_store,table_path', newLine: true, sort: true)

    emit:
    sdata = CREATE_SPATIALDATA.out.sdata // tuple(sample, zarr_store, table)
}
