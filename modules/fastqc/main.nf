// Performs quality check 
process FASTQC {

    label 'process_low_cpu'
    label 'process_low_memory'
    conda "${moduleDir}/env.yml"
    
    
    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.html"), emit: report
    tuple val(meta), path("*.zip"), emit: compressed_report
    path "versions.yml", emit: versions 

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def paired_end = meta.library_layout == 'PAIRED' ? true : false

    def thread_count = task.cpus
    def mem_mb = task.memory.toMega()
    args = "${args} --threads ${thread_count} --memory ${mem_mb} --outdir .".trim()
    
    def prefix = meta.prefix 

    def fastq_1 = null
    def fastq_2 = null
    // renaming files if necessary 
    if(paired_end){
        if(fastq instanceof List){
            if(fastq.size() == 2){
                fastq_1 = fastq[0]
                fastq_2 = fastq[1]
            }
            else{
                log.warn "Process FASTQC : Sample ${prefix} assigned FASTQ files >> ${fastq}."
                log.warn "Soft filtering with _1/_2 suffix was done to find the paired end file."
                fastq_1 = fastq.find { it =~ /_[1]\.(fq|fastq)(\.gz)?$/ }
                fastq_2 = fastq.find { it =~ /_[2]\.(fq|fastq)(\.gz)?$/ }
            }
        }
        else{
            error("${prefix} is paired-end but ${fastq} was given. Please consider verifying the metadata or exclude this sample")
        }
    }
    else{
        fastq_1 = fastq instanceof List ? fastq[0] : fastq
    }

    def already_prefixed_1 = fastq_1.name.startsWith(prefix)
    def extension1 = fastq_1.name.replaceAll(/^.*?\.(fastq|fq)/, '$1') 
    def link1 = already_prefixed_1 ? fastq_1.name : "${prefix}_raw_1.${extension1}"

    def already_prefixed_2 = false 
    def extension2 = ""
    def link2 = ""

    if(paired_end && fastq_2){
        already_prefixed_2 = fastq_2.name.startsWith(prefix)
        extension2 = fastq_2.name.replaceAll(/^.*?\.(fastq|fq)/, '$1')
        link2 = already_prefixed_2 ? fastq_2.name : "${prefix}_raw_2.${extension2}"
    }

    def fastqc_inputs = paired_end ? "${link1} ${link2}" : "${link1}"
    def stage = task.process.contains('POST') ? 'post_trim' : 'pre_trim'
    def target_dir_html = "${params.outdir}/${meta.parent}/fastqc/${stage}/report"
    def target_dir_zip  = "${params.outdir}/${meta.parent}/fastqc/${stage}/archive"

    """
    echo 'Starting Quality Check'
    echo 'Currently checking ${prefix}'

    shopt -s nullglob

    html_files=( ${target_dir_html}/${prefix}*.html )
    zip_files=( ${target_dir_zip}/${prefix}*.zip )

    if [ \${#html_files[@]} -gt 0 ] && [ \${#zip_files[@]} -gt 0 ]; then
        
        echo "[[FASTQC]] Files for ${prefix} already exist in publishDir. Copying to work dir..."
        
        cp ${target_dir_html}/${prefix}*.html ./
        cp ${target_dir_zip}/${prefix}*.zip ./

    else
        if [ "${already_prefixed_1}" = "false" ]; then
            ln -s ${fastq_1} ${link1}
        fi

        if [ "${paired_end}" = "true" ] && [ "${already_prefixed_2}" = "false" ]; then
            ln -s ${fastq_2} ${link2}
        fi

        fastqc \\
            ${fastqc_inputs} \\
            ${args} 

        echo 'Finished ${prefix}'
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastQC: \$(fastqc --version | sed 's/FastQC v//')
    END_VERSIONS
    """

}