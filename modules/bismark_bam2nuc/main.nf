// Bismark Nucleotide coverage report
process BISMARK_BAM2NUC {

    conda "${moduleDir}/env.yml"
    
    input:
    tuple val(meta), path(bam)
    path genome

    output:
    tuple val(meta), path("*.nucleotide_stats.txt"), emit: report
    path "versions.yml", emit: versions

    // look at ext.when to launch of not this process
    when:
    task.ext.when == null || task.ext.when

    script: 

    // get additional args if specified
    def args = task.ext.args ?: ''
    def prefix = meta.prefix
    args = "${args} --genome_folder ${genome} ${bam}"
    def target_dir_report = "${params.outdir}/${meta.parent}/bismark/nucleotide_coverage/report"

    """
    echo 'Starting Nucleotide coverage report ${prefix}'

    shopt -s nullglob
    report_files=( ${target_dir_report}/${prefix}*.nucleotide_stats.txt )

    if [ \${#report_files[@]} -gt 0 ]; then
        
        echo "[[BISMARK_BAM2NUC]] File for ${prefix} already exists in publishDir. Copying to work dir..."
        
        cp ${target_dir_report}/${prefix}*.nucleotide_stats.txt ./
        
    else
    
        bam2nuc \\
            ${args}
    fi

    echo 'Finished Nucleotide coverage report ${prefix}'


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bam2nuc: \$(bam2nuc --version | sed 's/.*version v//' | sed 's/ .*//')
    END_VERSIONS
    """

}