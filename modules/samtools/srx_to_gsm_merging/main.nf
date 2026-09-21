process SRX_TO_GSM_MERGING {

    label 'process_low_cpu'
    label 'process_low_memory'
    conda "${moduleDir}/env.yml"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.merged.bam"), emit: bam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    
    def single_run_check = meta.library_layout == 'PAIRED' ? (params.run_single_end ? true : false) : false    
    def target_dir_bam = "${params.outdir}/${meta.parent}/SRX_GSM_merged/bam"
    def prefix = meta.prefix 
    
    if(single_run_check){

        def bam_1 = bam.findAll { it =~ /_1.*\.bam$/ }
        def bam_2 = bam.findAll { it =~ /_2.*\.bam$/ }

        """
        shopt -s nullglob

        bam1=( ${target_dir_bam}/${prefix}_1.merged.bam )
        bam2=( ${target_dir_bam}/${prefix}_2.merged.bam )

        if [ \${#bam1[@]} -gt 0 ] && [ \${#bam2[@]} -gt 0 ]; then
            
            echo "[[SRX_TO_GSM_MERGING]] Files for ${prefix} already exist in publishDir. Copying to work dir..."
            
            cp ${target_dir_bam}/${prefix}_1.merged.bam .
            cp ${target_dir_bam}/${prefix}_2.merged.bam .
            
        else
            samtools merge -n -@ ${task.cpus - 1} ${prefix}_1.merged.bam ${bam_1}
            samtools merge -n -@ ${task.cpus - 1} ${prefix}_2.merged.bam ${bam_2}
        fi
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(samtools --version | head -1 | sed 's/samtools //')
        END_VERSIONS
        """
    }
    else{
        
        """
        shopt -s nullglob

        bam_files=( ${target_dir_bam}/${prefix}.merged.bam )

        if [ \${#bam_files[@]} -gt 0 ]; then
            
            echo "[[SRX_TO_GSM_MERGING]] File for ${prefix} already exists in publishDir. Copying to work dir..."
            
            cp ${target_dir_bam}/${prefix}.merged.bam ./
            
        else
            samtools merge -n -@ ${task.cpus - 1} ${prefix}.merged.bam ${bam}
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(samtools --version | head -1 | sed 's/samtools //')
        END_VERSIONS
        """
    }
}