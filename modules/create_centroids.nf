include { samplesFrom; publishedPath } from './samplesheet'

def centroidStem(sample) {
    params.group_by ? "${sample}.${params.group_by}.centroids" : "${sample}.centroids"
}

process CREATE_CENTROIDS {
    tag "${sample}"

    publishDir { "${params.outdir}/${sample}/create_centroids" }, mode: 'copy'

    // The store is carried, not staged: this step reads only the table.
    input:
    tuple val(sample), val(zarr_store), path(table)
    path 'timer.py'

    output:
    tuple val(sample), path("${centroidStem(sample)}.h5ad"), emit: centroids
    path "${centroidStem(sample)}.timing.tsv", emit: timings

    script:
    def centroidArgs = ["--sample ${sample}", "--table_path ${table}", "--outdir ."]
    if (params.group_by) centroidArgs << "--group_by ${params.group_by}"
    """
    create_centroids.py ${centroidArgs.join(' ')}
    """

    stub:
    """
    touch ${centroidStem(sample)}.h5ad
    touch ${centroidStem(sample)}.timing.tsv
    """
}

workflow create_centroids {
    take:
    // a samplesheet, or tuple(sample, zarr_store, table) per sample
    input

    main:
    def ch_tables = samplesFrom(input, ['sample', 'zarr_store', 'table_path'])

    CREATE_CENTROIDS(ch_tables, file("${projectDir}/bin/timer.py"))

    def sheet = params.group_by ? "create_centroids_${params.group_by}" : 'create_centroids'

    // Forwards the store and table it read beside the centroids it wrote: a report wants all
    // three, and only this step knows which table the centroids came from.
    CREATE_CENTROIDS.out.centroids
        .join(ch_tables)
        .map { sample, centroids, zarr_store, table ->
            "${sample},${publishedPath(input, sample, zarr_store)}," +
            "${publishedPath(input, sample, table)}," +
            "${params.outdir}/${sample}/create_centroids/${centroids.name}"
        }
        .collectFile(name: "${sheet}_samplesheet.csv", storeDir: params.outdir,
                     seed: 'sample,zarr_store,table_path,centroid_path', newLine: true, sort: true)

    emit:
    centroids = CREATE_CENTROIDS.out.centroids
}
