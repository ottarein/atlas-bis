// Aggregate quality check report  
process MULTIQC {

    conda "${moduleDir}/env.yml"
    
    input:
    tuple val(id), path(reports) 
    path multiqc_config

    output:
    path '*.html', emit: report
    path '*_data', emit : data
    path "versions.yml", emit: versions

    when: 
    task.ext.when == null || task.ext.when

    script: 
    """
    echo 'Initiate Multiqc : aggregation of QC report'
        
    multiqc . \\
        --config ${multiqc_config}

    echo 'Finished Multiqc'

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            multiqc: \$(multiqc --version | sed -e "s/multiqc, version //g")
    END_VERSIONS
    """
}