process TOOLS_VERSIONS {

    input:
    path('versions/*')

    output:
    path "software_versions.yml"

    script:
    """
    # Sort and deduplicate, keeping one entry per tool
    cat versions/* | sort -u > software_versions.yml
    """
}