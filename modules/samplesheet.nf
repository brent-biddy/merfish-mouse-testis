// With a column list, sample first then one file per remaining column, in the order given.
// Without one, sample and every path its row named, which is what a notebook needs: it
// declares what it wants by globbing, so render cannot name the columns.
def samplesFrom(input, List columns = null) {
    // steps.nf passes a samplesheet. Quotes honoured, so a value may hold a comma -- step 1a's
    // label_offset is ROWS,COLUMNS -- which splitCsv otherwise splits on.
    if (input instanceof Path || input instanceof String) {
        def ch_samples = channel.fromPath(input)
            .splitCsv(header: true, quote: '"')
            .map { row ->
                def cols = columns ?: ['sample'] + (row.keySet() - 'sample').toList()
                cols.each { c -> if (!row[c]) error "Samplesheet row missing '${c}': ${row}" }
                def paths = cols.tail().collect { file(row[it]) }
                columns ? [row.sample] + paths : tuple(row.sample, paths)
            }

        return ch_samples
    }
    // main.nf passes a channel already carrying one tuple per sample
    else {
        return input
    }
}

// Where a step's input was published, for the handoff sheet. A samplesheet's path already is
// that. A chained run's is a work dir the run deletes, so rebuild it: modules publish to
// ${params.outdir}/<sample>/<step>, and <step> is the middle field of <sample>.<step>.<ext>.
def publishedPath(input, sample, path) {
    if (input instanceof Path || input instanceof String) {
        return path.toString()
    }
    def step = path.name.substring("${sample}".length() + 1, path.name.lastIndexOf('.'))
    return "${params.outdir}/${sample}/${step}/${path.name}"
}
