process FASTQDL {

    label 'process_single'
    conda "${moduleDir}/env.yml"
        
    input:
    val(meta)
    val(sub_samplesheet) 

    output:
    tuple val(meta), path("*.{fq,fastq}{,.gz}"), emit : fastq
    path "versions.yml", emit: versions 

    when:
    (task.ext.when == null || task.ext.when) && \
    (sub_samplesheet == null || sub_samplesheet == [] || meta.run_accession in sub_samplesheet)

    script:
    def accession = meta.run_accession
    def target_dir = "${params.outdir}/fastqdl/fastq/"
    
    """
    shopt -s nullglob

    fastq_files=( ${target_dir}/${accession}*.{fastq,fq}{.gz,} )

    if [ \${#fastq_files[@]} -gt 0 ]; then
        
        echo "[[FASTQDL]] Files for ${accession} already exist in publishDir. Copying to work dir..."
        
        cp ${target_dir}/${accession}*.{fastq,fq}{.gz,} ./
        
    else
        echo "Downloading ${accession}"
        fastq-dl --accession ${accession}
        echo "Finished ${accession}"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastq-dl: \$(fastq-dl --version | sed 's/FastqDL v//')
    END_VERSIONS
    """
}