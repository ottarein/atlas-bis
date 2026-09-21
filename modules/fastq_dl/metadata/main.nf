// Retrieve all input metadata
process METADATA {

    label 'process_medium_cpu' 
    label 'process_lowX_memory' 
    conda "${moduleDir}/env.yml"
    
    input:
    val(samplesheet)

    output:
    path "*.tsv", emit: metadata
    path "versions.yml", emit: versions 

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = "--only-download-metadata"
    def accession = samplesheet.accession
    """
    FASTQDL_REGEX="^([EDS]R[PSXR][0-9]+|PRJ(NA|EB|DB)[0-9]+|SAM[DEN][A-Z0-9]+)\$"

    echo "Checking the format of the accession :"
    if [[ "${accession}" =~ \$FASTQDL_REGEX ]]; then 

        echo "The accession foormat is natively supported by fastql-dl; continue"
        echo "Retrieving metadata for ${accession}"

        fastq-dl \\
            --accession ${accession} \\
            ${args}

    elif [[ "${accession}" =~ ^GSE[0-9]+\$ ]]; then 

        echo "The accession format is detected as GSE from GEO, it will be converted :"
        SRP=\$(pysradb gse-to-srp "${accession}" | awk 'NR>1 {print \$2}' | head -n 1)
        echo "Convertion of ${accession} to \${SRP}"

        echo "Retrieving metadata for \${SRP}"

        fastq-dl \\
            --accession \${SRP} \\
            ${args}


    elif [[ "${accession}" =~ ^GSM[0-9]+\$ ]]; then

        echo "The accession format is detected as GSM from GEO, it will be converted :"
        SRR=\$(pysradb gsm-to-srr "${accession}" | awk 'NR>1 {print \$2}' | head -n 1)
        echo "Convertion of ${accession} to \${SRR}"

        echo "Retrieving metadata for \${SRR}"

        fastq-dl \\
            --accession \${SRR} \\
            ${args}
        
    else
        echo "Unsupported or unrecognized accession format: ${accession}" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastq-dl: \$(fastq-dl --version | sed 's/FastqDL v//')
    END_VERSIONS

    """
}