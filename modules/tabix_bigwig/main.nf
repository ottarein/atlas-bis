// Tabix indexing
process TABIX_BIGWIG {

    label 'process_medium_cpu'
    label 'process_lowX_memory'
    conda "${moduleDir}/env.yml"
    
    input:
    tuple val(meta), path(bedgraph)
    path genome

    output:
    tuple val(meta), path("*.sorted.bedGraph.gz"), emit: sorted_bedgraph
    tuple val(meta), path("*.bedGraph.gz.tbi"), emit: index
    tuple val(meta), path("*.bw*"), emit: bigwig
    path "versions.yml", emit: versions

    // look at ext.when to launch of not this process
    when:
    task.ext.when == null || task.ext.when

    script: 
    def dir_bedgraph_tabix = "${params.outdir}/${meta.parent}/bedgraph_tabix"
    def dir_bigwig = "${params.outdir}/${meta.parent}/bigwig"
    def prefix = meta.prefix
    def suffix = bedgraph.name.contains("UCSC") ? ".UCSC" : ""
    def name = "${prefix}${suffix}"

    """
    echo 'Starting Tabix indexing and BigWig generation for ${name}'

    shopt -s nullglob

    bg_files=( ${dir_bedgraph_tabix}/${name}.sorted.bedGraph.gz )
    tbi_files=( ${dir_bedgraph_tabix}/${name}.sorted.bedGraph.gz.tbi )
    bw_files=( ${dir_bigwig}/${name}.bw )

    if [ \${#bg_files[@]} -gt 0 ] && [ \${#tbi_files[@]} -gt 0 ] && [ \${#bw_files[@]} -gt 0 ]; then
       
        echo "[[TABIX_BIGWIG]] Files for ${name} already exist in publishDir. Copying to work dir..."
        
        cp ${dir_bedgraph_tabix}/${name}.sorted.bedGraph.gz ./
        cp ${dir_bedgraph_tabix}/${name}.sorted.bedGraph.gz.tbi ./
        cp ${dir_bigwig}/${name}.bw .

    else
    
        zcat ${bedgraph} | \\
            grep -v '^track' | \\
            LC_ALL=C sort -k1,1V -k2,2n -S ${task.memory.giga}G --parallel=${task.cpus} -T ./ > ${name}.sorted.bedGraph

        bedGraphToBigWig ${name}.sorted.bedGraph ${genome}/genome.chrom.sizes ${name}.bw

        bgzip ${name}.sorted.bedGraph

        tabix -p bed ${name}.sorted.bedGraph.gz
    fi

    echo 'Finished Tabix indexing and BigWig generation for ${name}'

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sort: \$(sort --version | head -n 1 | sed 's/sort (GNU coreutils) //')
            bgzip: \$(bgzip --version | head -n 1 | sed 's/htslib (bgzip) //')
            tabix: \$(tabix --version | head -n 1 | sed 's/htslib (tabix) //')
            ucsc-bedgraphtobigwig: \$(bedGraphToBigWig 2>&1 | head -n 1 | awk '{print \$2}')
    END_VERSIONS
    """

}