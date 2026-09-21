/*
 * Include modules
 */
include { METADATA } from "./modules/fastq_dl/metadata/main"
include { FASTQDL } from "./modules/fastq_dl/fastq/main"
include { FASTQC as PRE_TRIM_FASTQC } from "./modules/fastqc/main"
include { FASTQC as POST_TRIM_FASTQC } from "./modules/fastqc/main"
include { MULTIQC as MULTIQC_ACCESSION } from "./modules/multiqc/main"
include { MULTIQC as MULTIQC_STAGE } from "./modules/multiqc/main"
include { MULTIQC as MULTIQC_POST } from "./modules/multiqc/main"
include { TRIMGALORE } from "./modules/trimgalore/main"
include { BISMARK_GENOME_PREPARATION } from "./modules/bismark_genome_preparation/main"
include { BISMARK_ALIGNMENT } from "./modules/bismark_alignment/main"
include { BISMARK_DEDUPLICATION } from "./modules/bismark_deduplication/main"
include { BISMARK_METHYLATION_EXTRACTION } from "./modules/bismark_methylation_extraction/main"
include { BISMARK_BAM2NUC } from "./modules/bismark_bam2nuc/main"
include { SRX_TO_GSM_MERGING } from "./modules/samtools/srx_to_gsm_merging/main"
include { PAIRED_BAM_MERGING } from "./modules/samtools/paired_bam_merging/main"
include { QUALIMAP } from "./modules/qualimap/main"
include { TABIX_BIGWIG } from "./modules/tabix_bigwig/main"


// Create a function that foes into the metadata in be very beginning and search for any GSM+d so that it extract and create a new column called sample_GSM
// You need to assume that in the metadata there should be a GSM value and if not then you cant do SRR GSM merging and also you assume that if it exist then it is 
// of the sample and not some other information

// Change the publishDir for alignment to be StoreDir
// Change the genome_preparation publish dir to be storedir
// Change the methylation extraction publish_dir to be store_dir


/*
 * Pipeline
 */
workflow {

    // Dictionnary of accession type of what you should found in metadata
    def prefix_to_columns = [
    'SRR'   : "run_accession",
    'SRX'   : "experiment_accession",
    'GSM'   : "library_name;experiment_alias;sample_alias",
    'SAMN'  : "sample_accession", 
    'SRS'   : "secondary_sample_accession",
    'GSE'   : "study_alias",
    'PRJNA' : "study_accession",
    'SRP'   : "secondary_study_accession",
    'SRA'   : "submission_accession"
    ]

    // input samplesheet with all accession numbers
    ch_samplesheet = channel.fromPath(params.samplesheet, checkIfExists: true)
    .splitCsv(header : true, sep : ',')

    
    // Channel to collect all tools versions into a single channel
    def ch_versions = channel.empty()
    // Channel to collect all report into a single channel
    def ch_report = channel.empty()

    // Retrieve METADATA
    // if a metadata file is provided by the user, skip retrieval of metadata entirely
    // this is useful when raw metadata needs some reformatting before being use
    // it wasn't tested in the case of multiple larges studies, use with care
    if (params.metadata.path) {
        ch_metadata = channel.fromPath(params.metadata.path)
            .splitCsv(header: true, sep: '\t')
    } 
    else {
        METADATA(ch_samplesheet)
        ch_versions = ch_versions.mix(METADATA.out.versions.first())

        ch_metadata = METADATA.out.metadata
            .collectFile(name: 'combined_metadata.tsv', keepHeader: true)
            .splitCsv(header: true, sep: '\t')
    }

    // When in initialization, e.g simply retrieve the metadata, everything in the pipeline is skipped
    // If not, the pipeline is launched hereafter is used, the metadata is always retrieved being in Initialization
    // is only to get a feeling of what the samples are and complete correctly the samplesheet to give to the pipeline
    if (!params.init) {

        // define a subset of accession number from the metadata corresponding to the samples to download and study
        // do not confuse with "being a subset of the pipeline samplesheet" which is given to the pipeline
        // here the samplesheet refers to the metadata, so a subset of accession number inside the metadata
        // this files is simply a list of accession number in each line
        def sub_samplesheet = []
        if (params.sub_samplesheet) {
            sub_samplesheet = file(params.sub_samplesheet)
            .readLines()
            .collect { line -> line.trim().replaceAll('\\r', '') }
            .findAll { line -> line }
        }

        // retrieving and downloading the fastq files
        FASTQDL(ch_metadata, sub_samplesheet)
        ch_versions = ch_versions.mix(FASTQDL.out.versions.first())

        // look at the samplesheet given to the pipeline, 
        // infers the correct manipulation of the metadata in the form of a map and 
        // extract useful information for the pipeline to use after 
        def ch_rules_map = ch_samplesheet
        .reduce([:]) { map, row ->
            def acc = row.accession
            if (!map.containsKey(acc)) {
                map[acc] = [:]
                def matcher = acc =~ /^([A-Za-z]+)/
                def targetCols = []
                if (matcher.find()) {
                    def prefix = matcher.group(1).toUpperCase()
                    targetCols = prefix_to_columns[prefix] ?: []
                }
                
                map[acc]['type'] = targetCols
                map[acc]['stage'] = row.stage

            }
            
            def customsList = row.customs ? row.customs.split('\\$;\\$', -1)*.trim() : []
            def extractList = row.extract ? row.extract.split('\\$;\\$', -1)*.trim() : []
            def replaceList = row.replace ? row.replace.split('\\$;\\$', -1)*.trim() : []
            
            customsList.eachWithIndex { customCol, i ->
                if (customCol) { 
                    def ext = (i < extractList.size()) ? extractList[i] : ""
                    def rep = (i < replaceList.size()) ? replaceList[i] : ""
                    
                    map[acc][customCol] = [ extract: ext, replace: rep ]
                }
            }
            
            return map
        }

        def ch_metadata_prefix = FASTQDL.out.fastq
        .combine(ch_rules_map)
        .map { meta, fastq, rule_map -> 
            def final_prefix = ""
            def parent = ""
            def gsm = ""
            def skip_srx = false

            def gsm_cols = prefix_to_columns['GSM'].split(';')
            def gsm_values = gsm_cols.collect { col_name -> 
                meta[col_name.trim()] 
            }
            def valid_gsm = gsm_values.findAll { it != null && it.toString().trim() != '' }
            def unique_gsms = valid_gsm.unique()
            if (valid_gsm.isEmpty()) {
                skip_srx = true
                log.warn "Validation warning: No GSM values found inside the metadata for ${meta.accession}. SRX_TO_GSM will be altogether skipped."
            }
            else{
                if (unique_gsms.size() > 1) {
                log.warn "Validation warning: Conflicting GSM values found: ${unique_gsms}. Defaulting to the first one: ${unique_gsms[0]}"
                }
                gsm = unique_gsms[0]
            }         

            for(def entry in rule_map){
                
                def rule_acc = entry.key
                def rule_data = entry.value
                def type_col = rule_data.type

                if(meta[type_col]==rule_acc){
                    
                    def temp_list = []
                    rule_data.each { col_name, col_rules -> 
                        if (col_name != 'type' && col_name != 'stage') {
                            
                            def val = meta[col_name]
                            if(val!=null){
                                val = val.toString()
                                if(col_rules.extract){

                                    def cleaned_extract = col_rules.extract.replaceAll('^"|"$', '')
                                    def matcher = val =~ cleaned_extract
                                    if(matcher.find() && matcher.groupCount() >= 1){

                                        val = matcher.group(1)
                                    }
                                }
                                if(col_rules.replace){

                                    def parts = col_rules.replace.split(java.util.regex.Pattern.quote('$:$'))

                                    if(parts.size()==2){
                                        
                                        def origin = parts[0].trim().replaceAll('"', '')
                                        def target = parts[1].trim().replaceAll('"', '')

                                        val = val.replaceAll(origin, target)

                                    }
                                }

                                temp_list << val
                            }
                        }
                    }
                
                final_prefix = temp_list.join("_").replaceAll(' ', '_')
                parent = rule_acc
                break
                }
            }

            def new_meta = meta + [prefix : final_prefix, parent : parent, gsm : gsm, skip_srx : skip_srx]

            return [new_meta, fastq]
        }

        PRE_TRIM_FASTQC(ch_metadata_prefix)
        ch_versions = ch_versions.mix(PRE_TRIM_FASTQC.out.versions.first())
        ch_report = ch_report.mix(PRE_TRIM_FASTQC.out.compressed_report)

        TRIMGALORE(ch_metadata_prefix)
        ch_versions = ch_versions.mix(TRIMGALORE.out.versions.first())
        ch_report = ch_report.mix(TRIMGALORE.out.report)

        ch_metadata_prefix = TRIMGALORE.out.fastq
        POST_TRIM_FASTQC(ch_metadata_prefix)
        ch_versions = ch_versions.mix(POST_TRIM_FASTQC.out.versions.first())
        ch_report = ch_report.mix(POST_TRIM_FASTQC.out.compressed_report)

        multiqc_config = file("${projectDir}/modules/multiqc/multiqc_config.yml")

        if(!params.multiqc.break_parent){
            def ch_report_per_study = ch_report
            .map{ meta, report -> 
                def study = [ id : meta.parent ]
                return tuple(study, report)    
            }
            .groupTuple()
            .map{id, rep -> 
                def report = rep.flatten().collect()
                return tuple(id, report)
            }
            MULTIQC_ACCESSION(ch_report_per_study, multiqc_config)
            ch_versions = ch_versions.mix(MULTIQC_ACCESSION.out.versions)
        }
        else{
            MULTIQC_ACCESSION(ch_report, multiqc_config)
            ch_versions = ch_versions.mix(MULTIQC_ACCESSION.out.versions)
        }
        
        if(params.multiqc.per_stage){
            def ch_per_stage_report = ch_report
            .combine(ch_rules_map)
            .map{ meta, report, rule_map -> 
            def id = rule_map[meta.parent].stage
            def meta_stage = "${meta.parent}/${meta[id].toLowerCase()}"
            def study = [ id : meta_stage ]
            return tuple(study, report)    
            }
            .groupTuple()
            .map{id, rep -> 
            def report = rep.flatten().collect()
            return tuple(id, report)
            }
                
            MULTIQC_STAGE(ch_per_stage_report, multiqc_config)
            ch_versions = ch_versions.mix(MULTIQC_STAGE.out.versions)
        }

        // Bismark Pipeline
        def genome = channel.fromPath(params.genome).first()
        def ch_genome_index
        
        if (params.bismark_genome_preparation.skip) {
            ch_genome_index = genome
        } else {
            BISMARK_GENOME_PREPARATION(genome)
            ch_versions = ch_versions.mix(BISMARK_GENOME_PREPARATION.out.versions)
            ch_genome_index = BISMARK_GENOME_PREPARATION.out.output_dir
        }

        // Alignement
        BISMARK_ALIGNMENT(TRIMGALORE.out.fastq, ch_genome_index)
        ch_versions = ch_versions.mix(BISMARK_ALIGNMENT.out.versions.first())
        ch_report = ch_report.mix(BISMARK_ALIGNMENT.out.report)

        def ch_srx_to_gsm_skip = BISMARK_ALIGNMENT.out.bam
        .map { meta, bam -> meta.skip_srx }
        .collect()
        .map { flags -> flags.contains(true) }

        ch_alignment_evaluated = BISMARK_ALIGNMENT.out.bam
        .combine(ch_srx_to_gsm_skip)
        .branch { meta, bam, global_skip ->
            // If global_skip is true, bypass grouping for everyone.
            // We return just [meta, bam] so the boolean is removed.
            skip_merge: global_skip == true
                return [meta, bam]
            // If false, send everyone to be grouped.
            do_merge: global_skip == false
                return [meta, bam]
        }
    
        ch_grouped_bam = ch_alignment_evaluated.do_merge
        .map { meta, bam -> [ meta.gsm, meta, bam ] }// Grouping by GSM 
        .groupTuple()
        .map { gsm, metas, bams ->
            def group_meta = metas.first() + [id: gsm]
            [ group_meta, bams.flatten() ]
        }

        ch_grouped_bam
        .branch { meta, bams ->
            def count = bams instanceof List ? bams.size() : 1
            def min_to_merge = (params.run_single_end && meta.library_layout == 'PAIRED') ? 3 : 2
            
            merge: count >= min_to_merge
            skip:  true
                return [meta, (meta.library_layout == 'SINGLE' && bams instanceof List) ? bams[0] : bams]
        }
        .set { ch_srx_gsm_branched }

        SRX_TO_GSM_MERGING(ch_srx_gsm_branched.merge)
        ch_versions = ch_versions.mix(SRX_TO_GSM_MERGING.out.versions.first())

        ch_fusion_srx_gsm = ch_alignment_evaluated.skip_merge
        .mix(ch_srx_gsm_branched.skip)
        .mix(SRX_TO_GSM_MERGING.out.bam)

        def ch_qualimap = ch_fusion_srx_gsm.transpose()

        QUALIMAP(ch_qualimap)
        ch_report = ch_report.mix(QUALIMAP.out.report)
        ch_versions = ch_versions.mix(QUALIMAP.out.versions.first())

        BISMARK_DEDUPLICATION(ch_fusion_srx_gsm)
        ch_report = ch_report.mix(BISMARK_DEDUPLICATION.out.report)
        ch_versions = ch_versions.mix(BISMARK_DEDUPLICATION.out.versions.first())

        if(params.run_single_end){
            BISMARK_DEDUPLICATION.out.bam
            .branch { meta, bams ->
                merge: meta.library_layout == 'PAIRED' && bams instanceof List && bams.size() > 1
                skip:  true
                    return [meta, bams instanceof List ? bams[0] : bams]
            }
            .set { ch_paired_to_single_branched }
            PAIRED_BAM_MERGING(ch_paired_to_single_branched.merge)
            ch_versions = ch_versions.mix(PAIRED_BAM_MERGING.out.versions.first())

            ch_deduplicated_bam = PAIRED_BAM_MERGING.out.bam.mix(ch_paired_to_single_branched.skip)
        }
        else{
            ch_deduplicated_bam = BISMARK_DEDUPLICATION.out.bam
        }
        
        // Methylation Extraction
        BISMARK_METHYLATION_EXTRACTION(ch_deduplicated_bam, genome)
        ch_report = ch_report.mix(BISMARK_METHYLATION_EXTRACTION.out.report)
        ch_report = ch_report.mix(BISMARK_METHYLATION_EXTRACTION.out.mbias)
        ch_report = ch_report.mix(BISMARK_METHYLATION_EXTRACTION.out.cytosine_context_summary)
        ch_versions = ch_versions.mix(BISMARK_METHYLATION_EXTRACTION.out.versions.first())

        //Nucleotide coverage report
        BISMARK_BAM2NUC(ch_deduplicated_bam, genome) 
        ch_report = ch_report.mix(BISMARK_BAM2NUC.out.report)
        ch_versions = ch_versions.mix(BISMARK_BAM2NUC.out.versions.first())  

        TABIX_BIGWIG(BISMARK_METHYLATION_EXTRACTION.out.bedgraph.transpose(), genome)   
        ch_versions = ch_versions.mix(TABIX_BIGWIG.out.versions.first())    

        if(!params.multiqc.break_parent){
            def ch_report_per_study_post = ch_report
            .map{ meta, report -> 
                def study = [ id : meta.parent ]
                return tuple(study, report)    
            }
            .groupTuple()
            .map{id, rep -> 
                def report = rep.flatten().collect()
                return tuple(id, report)
            }
            MULTIQC_POST(ch_report_per_study_post, multiqc_config)
            ch_versions = ch_versions.mix(MULTIQC_POST.out.versions)
        }
        else{
            MULTIQC_POST(ch_report, multiqc_config)
            ch_versions = ch_versions.mix(MULTIQC_POST.out.versions)
        }

        // Saving tools versions
        ch_versions
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info")
    }
}