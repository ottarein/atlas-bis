// Bismark genome preparation
process BISMARK_GENOME_PREPARATION {

    label 'process_high_cpu'
    label 'process_medium_memory'
    conda "${moduleDir}/env.yml"

    input:
    // It is important to know that this is a symlink
    path(genome, stageAs: "BismarkIndex")

    output:
    path "BismarkIndex", emit: output_dir
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    echo "Checking if Bismark Index already exists..."
    
    if [ -d "${genome}/Bisulfite_Genome" ] && [ -f "${genome}/genome.chrom.sizes" ]; then
        
        echo "[[BISMARK_GENOME_PREPARATION]] Directory and file for genome already exist in the Genome folder. Proceeding..."
        
    else
        echo "Starting Bismark Genome Preparation..."
        
        bismark_genome_preparation \\
            --verbose \\
            --parallel ${task.cpus} \\
            BismarkIndex

        FASTA=\$(find BismarkIndex/ -maxdepth 1 -name '*.fa' -o -name '*.fasta' -o -name '*.fa.gz'| head -1)

        samtools faidx \${FASTA}

        cut -f1,2 \${FASTA}.fai > BismarkIndex/genome.chrom.sizes
        
        echo "Finished Bismark Genome Preparation"
    fi

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bismark: \$(bismark --version | head -n 1 | sed 's/.*Bismark Version: v//')
            samtools: \$(samtools --version | head -n 1 | sed 's/samtools //')
    END_VERSIONS
    """
}