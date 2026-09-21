process BISMARK_DEDUPLICATION {

    label 'process_low_cpu'
    label 'process_low_memory'
    conda "${moduleDir}/env.yml"
    
    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.deduplicated.bam"), emit: bam
    tuple val(meta), path("*.deduplication_report.txt"), emit: report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when
      
    script: 
    def args = task.ext.args ?: ''
    def prefix = meta.prefix
    def run_type = meta.library_layout == 'PAIRED' ? '--paired' : '--single'
    def target_dir_bam = "${params.outdir}/${meta.parent}/bismark/deduplication/bam"
    def target_dir_report = "${params.outdir}/${meta.parent}/bismark/deduplication/report"

    def bismark_cmd = ""
    def bypass_check = ""

    if(meta.library_layout == 'PAIRED' && params.run_single_end){
        
        bypass_check = """bam1=( ${target_dir_bam}/${prefix}*{_1,val_1}*.deduplicated.bam )
        bam2=( ${target_dir_bam}/${prefix}*{_2,val_2}*.deduplicated.bam )
        txt1=( ${target_dir_report}/${prefix}*{_1,val_1}*.deduplication_report.txt )
        txt2=( ${target_dir_report}/${prefix}*{_2,val_2}*.deduplication_report.txt )

        if [ \${#bam1[@]} -gt 0 ] && [ \${#bam2[@]} -gt 0 ] && [ \${#txt1[@]} -gt 0 ] && [ \${#txt2[@]} -gt 0 ]; then
        """

        bismark_cmd = """for bam_file in ${bam}; do
            echo "Processing \${bam_file}..."
            deduplicate_bismark \\
                ${args} \\
                --single \\
                \${bam_file}
        done
        """
    }
    else {

        bypass_check = """bam_files=( ${target_dir_bam}/${prefix}*.deduplicated.bam )
        txt_files=( ${target_dir_report}/${prefix}*.deduplication_report.txt )

        if [ \${#bam_files[@]} -gt 0 ] && [ \${#txt_files[@]} -gt 0 ]; then
        """

        bismark_cmd = """deduplicate_bismark \\
            ${args} \\
            ${run_type} \\
            ${bam}
        """
    }

    """
    echo 'Starting deduplication ${prefix}'

    shopt -s nullglob

    ${bypass_check}
        
        echo "[[BISMARK_DEDUPLICATION]] Files for ${prefix} already exist in publishDir. Copying to work dir..."
        
        cp ${target_dir_bam}/${prefix}*.deduplicated.bam ./
        cp ${target_dir_report}/${prefix}*.deduplication_report.txt ./
        
    else
        ${bismark_cmd}
    fi

    echo 'Finished deduplication ${prefix}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deduplicate_bismark: \$(deduplicate_bismark --version | sed -n 's/.*Deduplicator Version: //p')
    END_VERSIONS
    """
}