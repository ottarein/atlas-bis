// Bismark alignment
process BISMARK_ALIGNMENT {

    label 'process_extended_cpu'
    label 'process_high_memory'
    conda "${moduleDir}/env.yml"
    
    input:
    tuple val(meta), path(fastq)
    path genome

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.txt"), emit: report
    path "versions.yml", emit: versions 

    // look at ext.when to launch of not this process
    when:
    task.ext.when == null || task.ext.when

    script: 

    // get additional args if specified
    def args = task.ext.args ?: ''
    // check paired-end sample
    def paired_end = meta.library_layout == 'PAIRED' ? true : false
    // check if the sample should be considered single-end either way
    def single_end_run = params.run_single_end ? true : false
    def prefix = meta.prefix

    // Construct the fastq arguments considering paired-end and single-end run arguments
    // if it is paired and need single-end run then the fastq will be handle afterwards thus an empty string here
    def fastq_ = paired_end ? (single_end_run ? "" : "-1 ${fastq[0]} -2 ${fastq[1]}") : fastq

    // Try to assign sensible bismark --multicore if not already set
    if(!args.contains('--multicore') && task.cpus){

        // Numbers based on recommendation by Felix for a typical mouse genome
        def ccore = 1
        def cpu_per_multicore = 3
        def mem_per_multicore = (13.GB).toBytes()
        if(args.contains('--non_directional')){
            cpu_per_multicore = 5
            mem_per_multicore = (18.GB).toBytes()
        }

        // How many multicore splits can we afford with the cpus we have?
        ccore = ((task.cpus as int) / cpu_per_multicore) as int

        // Check that we have enough memory
        try {
            def tmem = (task.memory as MemoryUnit).toBytes()
            def mcore = (tmem / mem_per_multicore) as int
            ccore = Math.min(ccore, mcore)
        } catch (Exception e) {
            log.warn "Error catched: ${e}"
            log.warn "Not able to define bismark align multicore based on available memory"
        }
        if(ccore > 1){
            args += " --multicore ${ccore}"
        }
    }
    
    def single_args = "${args} --genome ${genome}"
    def general_args = "${args} --genome ${genome} ${fastq_}"

    def target_dir_bam = "${params.outdir}/${meta.parent}/bismark/alignement/bam"
    def target_dir_report = "${params.outdir}/${meta.parent}/bismark/alignement/report"

    def bismark_cmd = ""
    def bypass_check = ""

    if(single_end_run && paired_end){
        
        bypass_check = """bam1=( ${target_dir_bam}/${prefix}*{_1,val_1}*.bam )
        bam2=( ${target_dir_bam}/${prefix}*{_2,val_2}*.bam )
        txt1=( ${target_dir_report}/${prefix}*{_1,val_1}*.txt )
        txt2=( ${target_dir_report}/${prefix}*{_2,val_2}*.txt )

        if [ \${#bam1[@]} -gt 0 ] && [ \${#bam2[@]} -gt 0 ] && [ \${#txt1[@]} -gt 0 ] && [ \${#txt2[@]} -gt 0 ]; then
        """

        bismark_cmd = """echo 'Paired-end sample will be processed individually following arguments : run_single_end'
        bismark \\
            ${single_args} \\
            ${fastq[0]}
        bismark \\
            ${single_args} \\
            ${fastq[1]}
        """

    } else {
        
        bypass_check = """bam_files=( ${target_dir_bam}/${prefix}*.bam )
        txt_files=( ${target_dir_report}/${prefix}*.txt )

        if [ \${#bam_files[@]} -gt 0 ] && [ \${#txt_files[@]} -gt 0 ]; then
        """

        bismark_cmd = """bismark \\
            ${general_args} 
        """
    }

    """
    echo 'Preparing alignment ${prefix}'

    shopt -s nullglob

    ${bypass_check}
       
        echo "[[BISMARK_ALIGNMENT]] Files for ${prefix} already exist in publishDir. Copying to work dir..."
        
        cp ${target_dir_bam}/${prefix}*.bam ./
        cp ${target_dir_report}/${prefix}*.txt ./
        
    else
        echo 'Starting Bismark alignment...'
        
        ${bismark_cmd}
        
        echo 'Finished alignment ${prefix}'
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Bismark: \$(bismark --version | sed 's/Bismark:Alignment v//' | head -n 1)
    END_VERSIONS
    """
}