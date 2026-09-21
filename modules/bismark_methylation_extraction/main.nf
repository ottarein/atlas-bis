process BISMARK_METHYLATION_EXTRACTION {

    label 'process_high_cpu'
    label 'process_medium_memory'
    conda "${moduleDir}/env.yml"
    
    input:
    tuple val(meta), path(bam)
    path genome

    output:
    tuple val(meta), path("*.bedGraph.gz"), emit: bedgraph
    tuple val(meta), path("*.txt.gz"), emit: methylation_calls
    tuple val(meta), path("*.cov.gz"), emit: coverage
    tuple val(meta), path("*_splitting_report.txt"), emit: report
    tuple val(meta), path("*.M-bias.txt"), emit: mbias
    tuple val(meta), path("*.cytosine_context_summary.txt"), emit: cytosine_context_summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: 

    def args = task.ext.args ?: ''
    def prefix = meta.prefix
    def paired_end = meta.library_layout == 'PAIRED' ? true : false 
    def paired = paired_end ? '--paired-end' : '--single-end'

    if(paired_end && params.run_single_end){
        paired = '--single-end'
    }

    if(!args.contains('--multicore') && task.cpus >= 6){
        args += " --multicore ${(task.cpus / 3) as int}"
    }
    if(!args.contains('--buffer_size') && task.memory?.giga > 6){
        args += " --buffer_size ${task.memory.giga - 2}G"
    }

    args = "${args} ${paired} --genome_folder ${genome} --bedGraph --ucsc --counts --gzip --report --cytosine_report ${bam}"

    def dir_bedgraph = "${params.outdir}/${meta.parent}/bismark/methylation_extraction/bedgraph"
    def dir_coverage = "${params.outdir}/${meta.parent}/bismark/methylation_extraction/coverage"
    def dir_report = "${params.outdir}/${meta.parent}/bismark/methylation_extraction/splitting_report"
    def dir_mbias = "${params.outdir}/${meta.parent}/bismark/methylation_extraction/mbias"
    def dir_calls = "${params.outdir}/${meta.parent}/bismark/methylation_extraction/methylation_calls"
    def dir_summary = "${params.outdir}/${meta.parent}/bismark/methylation_extraction/cytosine_context_summary"
    def expected_calls_count = params.bismark_alignment.args.contains('--non_directional') ? 13 : 7
    
    """
    echo 'Extracting Methylation ${prefix}'

    shopt -s nullglob

    calls_files=( ${dir_calls}/*${prefix}*.txt.gz )
    bedgraph_files=( ${dir_bedgraph}/*${prefix}*.bedGraph.gz )
    cov_files=( ${dir_coverage}/*${prefix}*.cov.gz )
    report_files=( ${dir_report}/*${prefix}*_splitting_report.txt )
    mbias_files=( ${dir_mbias}/*${prefix}*.M-bias.txt )
    summary_files=( ${dir_summary}/*${prefix}*cytosine_context_summary.txt )

    if [ \${#calls_files[@]} -eq ${expected_calls_count} ] && [ \${#bedgraph_files[@]} -gt 0 ] && [ \${#cov_files[@]} -gt 0 ] && [ \${#report_files[@]} -gt 0 ] && [ \${#mbias_files[@]} -gt 0 ] && [ \${#summary_files[@]} -gt 0 ]; then

        echo "[[BISMARK_METHYLATION_EXTRACTION]] Files for ${prefix} already exist in publishDir. Copying to work dir..."

        cp ${dir_bedgraph}/${prefix}*.bedGraph.gz ./
        cp ${dir_coverage}/${prefix}*.cov.gz ./
        cp ${dir_report}/${prefix}*_splitting_report.txt ./
        cp ${dir_mbias}/${prefix}*.M-bias.txt ./
        cp ${dir_calls}/${prefix}*.txt.gz ./
        cp ${dir_summary}/${prefix}*cytosine_context_summary.txt ./

    else
        bismark_methylation_extractor \\
            ${args} 
    fi

    echo 'Finished methylation calling ${prefix}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bismark_methylation_extractor: \$(bismark_methylation_extractor --version | sed 's/.*version v//' | sed 's/ .*//')
    END_VERSIONS
    """

}