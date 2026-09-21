// Performs alignment quality check 
process QUALIMAP {

    label 'process_medium_cpu'
    label 'process_medium_memory'
    conda "${moduleDir}/env.yml"
    
    
    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.sorted_stats"), emit: report
    path "versions.yml", emit: versions 

    when:
    task.ext.when == null || task.ext.when

    script:

    def mem_per_core_mb = (task.memory.mega / task.cpus).toInteger()

    def suffix = bam.name.replace(meta.prefix, "")
    def read_id = suffix.find(/val_[12]/) ?: suffix.find(/_[12]/) ?: ""
    def prefix = meta.library_layout == 'PAIRED' ? (params.run_single_end ? meta.prefix + read_id : meta.prefix) : meta.prefix
    def target_dir_qualimap = "${params.outdir}/${meta.parent}/bamqc_report"
    
    """
    echo "Starting quality check on ${prefix} alignment bam"

    shopt -s nullglob

    qualimap_files=( ${target_dir_qualimap}/${prefix}*.sorted_stats )

    if [ \${#qualimap_files[@]} -gt 0 ]; then
        
        echo "[[QUALIMAP]] Files for ${prefix} already exist in publishDir. Copying to work dir..."
        
        cp -r ${target_dir_qualimap}/${prefix}*.sorted_stats ./
        
    else
        SORT_ORDER=\$(samtools view -H "${bam}" | grep -m 1 "^@HD" | grep -o "SO:[a-zA-Z]*" | cut -d: -f2)

        if [ "\$SORT_ORDER" = "coordinate" ]; then
            echo "Bam file already coordinate sorted, proceeding with bamqc"
            
            qualimap bamqc \\
                -bam ${bam} \\
                -outdir ${prefix}_bamqc \\
                -nt ${task.cpus} \\
                -c 

        else
            samtools sort -@ ${task.cpus} -m ${mem_per_core_mb}M -o "${prefix}.sorted.bam" "${bam}"

            qualimap bamqc \\
                -bam ${prefix}.sorted.bam \\
                -outdir ${prefix}_bamqc \\
                -nt ${task.cpus} \\
                -c 
                
        fi
    fi

    echo "Finished quality check on ${prefix} alignment bam"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        qualimap: \$(qualimap --version 2>&1 | grep 'QualiMap v.' | sed 's/QualiMap v.//')
    END_VERSIONS
    """

}