process PAIRED_BAM_MERGING {

    label 'process_low_cpu'
    label 'process_low_memory'
    conda "${moduleDir}/env.yml"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.paired_merged.bam"), emit: bam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def target_dir_bam = "${params.outdir}/${meta.parent}/paired_merged/bam"
    def prefix = meta.prefix

    """
    shopt -s nullglob

    bam_files=( ${target_dir_bam}/${prefix}.paired_merged.bam )

    if [ \${#bam_files[@]} -gt 0 ]; then
       
        echo "[[PAIRED_BAM_MERGING]] File for ${prefix} already exists in publishDir. Copying to work dir..."
        
        cp ${target_dir_bam}/${prefix}.paired_merged.bam ./
        
    else
        echo "Merging Paired BAM files for ${prefix}..."
        
        samtools merge -n -@ ${task.cpus - 1} ${prefix}.paired_merged.bam ${bam}

        echo "Finished merging Paired BAM files for ${prefix}."
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}