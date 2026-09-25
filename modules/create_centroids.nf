include { samplesFrom } from './samplesheet'

def centroidStem(sample) {
    params.group_by ? "${sample}.${params.group_by}.centroids" : "${sample}.centroids"
}

process CREATE_CENTROIDS {
    tag "${sample}"

    publishDir { "${params.outdir}/${sample}/create_centroids" }, mode: 'copy'

    input:
    tuple val(sample), path(zarr)
    path 'timer.py'

    output:
    tuple val(sample), path("${centroidStem(sample)}.h5ad"), emit: centroids
    path "${centroidStem(sample)}.timing.tsv", emit: timings

    script:
    def centroidArgs = ["--sample ${sample}", "--path ${zarr}", "--outdir ."]
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
    // a samplesheet, or tuple(sample, zarr) per sample
    input

    main:
    def ch_zarrs = samplesFrom(input, ['sample', 'path'])

    // A samplesheet's path column is already the published location. A chained run's is the
    // work dir the run deletes, so rebuild it: modules publish to ${params.outdir}/<sample>/<step>,
    // and <step> is the middle field of <sample>.<step>.<ext>.
    def ch_sources = input instanceof Path || input instanceof String
        ? ch_zarrs.map { sample, zarr -> tuple(sample, zarr.toString()) }
        : ch_zarrs.map { sample, zarr ->
            def step = zarr.name.substring("${sample}".length() + 1, zarr.name.lastIndexOf('.'))
            tuple(sample, "${params.outdir}/${sample}/${step}/${zarr.name}")
        }

    CREATE_CENTROIDS(ch_zarrs, file("${projectDir}/bin/timer.py"))

    def sheet = params.group_by ? "create_centroids_${params.group_by}" : 'create_centroids'

    CREATE_CENTROIDS.out.centroids
        .join(ch_sources)
        .map { sample, centroids, source ->
            "${sample},${source},${params.outdir}/${sample}/create_centroids/${centroids.name}"
        }
        .collectFile(name: "${sheet}_samplesheet.csv", storeDir: params.outdir,
                     seed: 'sample,path,centroid_path', newLine: true, sort: true)

    emit:
    centroids = CREATE_CENTROIDS.out.centroids
}
