process TRIMGALORE {

    label 'process_medium_cpu'
    label 'process_low_memory'
    conda "${moduleDir}/env.yml"
    
    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*{3prime,5prime,trimmed,val}{,_1,_2}.fq.gz"), emit: fastq
    tuple val(meta), path("*.txt"), emit: report
    tuple val(meta), path("*.json"), emit: summary
    path "versions.yml", emit: versions 

    when:
    task.ext.when == null || task.ext.when

    script: 
    def args = task.ext.args ?: ''
    def paired_end = meta.library_layout == 'PAIRED' ? true : false

    def accession = meta.run_accession
    def prefix = meta.prefix 

    def fastq_1 = null
    def fastq_2 = null
    
    if(paired_end){
        if(fastq instanceof List){
            if(fastq.size() == 2){
                fastq_1 = fastq[0]
                fastq_2 = fastq[1]
            }
            else{
                log.warn "Process TRIMGALORE : Sample ${accession} assigned FASTQ files >> ${fastq}."
                log.warn "Soft filtering with _1/_2 suffix was done to find the paired end file."
                fastq_1 = fastq.find { it =~ /_[1]\.(fq|fastq)(\.gz)?$/ }
                fastq_2 = fastq.find { it =~ /_[2]\.(fq|fastq)(\.gz)?$/ }
            }
        }
        else{
            error("${accession} is paired-end but ${fastq} was given. Please consider verifying the metadata or exclude this sample")
        }
    }
    else{
        fastq_1 = fastq instanceof List ? fastq[0] : fastq
    }

    def extension1 = fastq_1.name.replaceAll(/^.*?\.(fastq|fq)/, '$1') 
    def link1 = "${prefix}_1.${extension1}"
    
    def extension2 = ""
    def link2 = ""

    if(paired_end && fastq_2){
        extension2 = fastq_2.name.replaceAll(/^.*?\.(fastq|fq)/, '$1')
        link2 = "${prefix}_2.${extension2}"
    }

    def fastq_input = paired_end ? "--paired ${link1} ${link2}" : "${link1}"

    def cores = 1
    if (task.cpus) {
        cores = (task.cpus as int) - 4
        if (!paired_end) {
            cores = (task.cpus as int) - 3
        }
        if (cores < 1) {
            cores = 1
        }
        if (cores > 8) {
            cores = 8
        }
    }

    args = "${args} --cores ${cores}" 
    def target_dir_fastq   = "${params.outdir}/${meta.parent}/trimgalore/trimmed_fastq"
    def target_dir_report  = "${params.outdir}/${meta.parent}/trimgalore/report"
    def target_dir_summary = "${params.outdir}/${meta.parent}/trimgalore/summary_command"

    """
    echo 'Begin Trimming'
    echo 'Currently trimming ${accession}'

    shopt -s nullglob

    fastq_files=( ${target_dir_fastq}/${prefix}*{3prime,5prime,trimmed,val}{,_1,_2}.fq.gz )
    report_files=( ${target_dir_report}/${prefix}*.txt )
    summary_files=( ${target_dir_summary}/${prefix}*.json )

    if [ \${#fastq_files[@]} -gt 0 ] && [ \${#report_files[@]} -gt 0 ] && [ \${#summary_files[@]} -gt 0 ]; then
       
        echo "[[TRIMGALORE]] Files for ${accession} already exist in publishDir. Copying to work dir..."
        
        cp ${target_dir_fastq}/${prefix}*.fq.gz ./
        cp ${target_dir_report}/${prefix}*.txt ./
        cp ${target_dir_summary}/${prefix}*.json ./
        
    else
        ln -s ${fastq_1} ${link1}
        if [ "${paired_end}" = "true" ]; then
            ln -s ${fastq_2} ${link2}
        fi

        trim_galore \\
            ${fastq_input} \\
            ${args} 

        echo 'Finished Trimming ${accession}'
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        TrimGalore!: \$(trim_galore --version | sed 's/TrimGalore! v//')
    END_VERSIONS
    """
}