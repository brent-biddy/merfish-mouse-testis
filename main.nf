#!/usr/bin/env nextflow

include { create_spatialdata      } from './modules/create_spatialdata'
include { cluster_spatialdata_gpu } from './modules/cluster_spatialdata_gpu'
include { annotate_celltypes      } from './modules/annotate_celltypes'
include { create_centroids        } from './modules/create_centroids'
include { render                  } from './modules/render'

workflow {
    create_spatialdata(file(params.samplesheet))
    cluster_spatialdata_gpu(create_spatialdata.out.sdata)
    annotate_celltypes(cluster_spatialdata_gpu.out.sdata)

    create_centroids(annotate_celltypes.out.sdata)

    create_centroids.out.centroids
        .join(annotate_celltypes.out.sdata)
        .map { sample, centroids, zarr_store, table -> tuple(sample, [zarr_store, table, centroids]) }
        .set { ch_report_samples } // tuple(sample, staged_paths)

    render(ch_report_samples, file("${projectDir}/notebooks/celltype_report.qmd"), params.to, false)
}
