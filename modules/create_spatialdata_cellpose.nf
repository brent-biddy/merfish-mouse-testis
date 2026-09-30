include { samplesFrom } from './samplesheet'

process CREATE_SPATIALDATA_CELLPOSE {
    tag "${sample}"

    publishDir { "${params.outdir}/${sample}/create_spatialdata_cellpose" }, mode: 'copy'

    input:
    tuple val(sample), path(region_dir), path(vpt_dir)
    path 'timer.py'

    output:
    tuple val(sample), path("${sample}.create_spatialdata_cellpose.zarr"),
          path("${sample}.create_spatialdata_cellpose.h5ad"), emit: sdata
    path "${sample}.create_spatialdata_cellpose.timing.tsv", emit: timings

    script:
    """
    create_spatialdata_cellpose.py --sample ${sample} --path ${region_dir} \\
        --vpt_path ${vpt_dir} --outdir .
    """

    stub:
    """
    mkdir -p ${sample}.create_spatialdata_cellpose.zarr
    touch ${sample}.create_spatialdata_cellpose.h5ad
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

    CREATE_SPATIALDATA_CELLPOSE.out.sdata
        .map { sample, zarr_store, table ->
            def published = "${params.outdir}/${sample}/create_spatialdata_cellpose"
            "${sample},${published}/${zarr_store.name},${published}/${table.name}"
        }
        .collectFile(name: 'create_spatialdata_cellpose_samplesheet.csv', storeDir: params.outdir,
                     seed: 'sample,zarr_store,table_path', newLine: true, sort: true)

    emit:
    sdata = CREATE_SPATIALDATA_CELLPOSE.out.sdata // tuple(sample, zarr_store, table)
}
