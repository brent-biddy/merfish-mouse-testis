include { samplesFrom; publishedPath } from './samplesheet'

process EXPORT_SPATIALDATA {
    tag "${sample}"

    publishDir { "${params.outdir}/${sample}/export_spatialdata" }, mode: 'copy'

    input:
    tuple val(sample), path(zarr_store), path(table)
    path 'timer.py'

    output:
    tuple val(sample), path("${sample}.export_spatialdata.zarr"), emit: zarr
    path "${sample}.export_spatialdata.timing.tsv", emit: timings

    script:
    """
    export_spatialdata.py --sample ${sample} --zarr_store ${zarr_store} \\
        --table_path ${table} --outdir .
    """

    stub:
    """
    mkdir -p ${sample}.export_spatialdata.zarr
    touch ${sample}.export_spatialdata.timing.tsv
    """
}

workflow export_spatialdata {
    take:
    // a samplesheet, or tuple(sample, zarr_store, table) per sample
    input

    main:
    def ch_sdata = samplesFrom(input, ['sample', 'zarr_store', 'table_path'])

    EXPORT_SPATIALDATA(ch_sdata, file("${projectDir}/bin/timer.py"))

    // The exported store holds the table too; table_path still names the one it came from,
    // so this sheet has the same shape as every other step's.
    EXPORT_SPATIALDATA.out.zarr
        .join(ch_sdata)
        .map { sample, exported, zarr_store, table ->
            "${sample},${params.outdir}/${sample}/export_spatialdata/${exported.name}," +
            "${publishedPath(input, sample, table)}"
        }
        .collectFile(name: 'export_spatialdata_samplesheet.csv', storeDir: params.outdir,
                     seed: 'sample,zarr_store,table_path', newLine: true, sort: true)

    emit:
    zarr = EXPORT_SPATIALDATA.out.zarr
}
