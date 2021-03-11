/*
  Helper processes
*/

output_docs = Channel.fromPath("$baseDir/docs/output.md")

process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r ${file(output_docs)} results_description.html
    """
}


process get_software_versions {
    output:
    path 'software_versions_mqc.yaml'

    script:
    """
    echo $workflow.manifest.version &> v_ngi_rnaseq.txt
    echo $workflow.nextflow.version &> v_nextflow.txt
    fastqc --version &> v_fastqc.txt
    cutadapt --version &> v_cutadapt.txt
    trim_galore --version &> v_trim_galore.txt
    STAR --version &> v_star.txt
    hisat2 --version &> v_hisat2.txt
    stringtie --version &> v_stringtie.txt
    preseq &> v_preseq.txt
    read_duplication.py --version &> v_rseqc.txt
    echo \$(bamCoverage --version 2>&1) > v_deeptools.txt
    featureCounts -v &> v_featurecounts.txt
    picard MarkDuplicates --version &> v_markduplicates.txt  || true
    samtools --version &> v_samtools.txt
    multiqc --version &> v_multiqc.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}


/*
  Helper functions
*/
def helpMessage() {
    log.info """
    =======================================================
                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\'
        |\\ | |__  __ /  ` /  \\ |__) |__         }  {
        | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                              `._,._,\'

     nf-core/rnaseq : RNA-Seq Best Practice v${workflow.manifest.version}
    =======================================================

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/rnaseq --reads '*_R{1,2}.fastq.gz' --genome GRCh37 -profile uppmax

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      --readPathFile                Tab-seperated file with sample names and path to the fastq files. (Used if --reads not provided.)
      -profile                      Configuration profile to use. uppmax / uppmax_modules / hebbe / docker / aws

    Additional quantification options:
      --skip_alignment              Skip alignment with STAR or HISAT2 (deafault: false)
      --run_salmon                  Runs transcript expression quantification (Salmon)
      --run_txrevise                Runs txrevise quantification (Salmon with custom reference transciptome)
      --run_leafcutter              Runs alternative splicing quantification  (LeafCutter)
      --run_exon_quant              Runs exon quantification (DEXseq)

    Options:
      --genome                      Name of iGenomes reference
      --singleEnd                   Specifies that the input is single end reads
    Strandedness:
      --forward_stranded            The library is forward stranded
      --reverse_stranded            The library is reverse stranded
      --unstranded                  The default behaviour

    References                      If not specified in the configuration file or you wish to overwrite any of the references.
      --star_index                  Path to STAR index
      --hisat2_index                Path to HiSAT2 index
      --fasta                       Path to Fasta reference
      --gtf                         Path to GTF file
      --gff                         Path to GFF3 file
      --bed12                       Path to bed12 file
      --saveReference               Save the generated reference files the the Results directory.
      --saveTrimmed                 Save trimmed FastQ file intermediates
      --saveAlignedIntermediates    Save the BAM files from the Aligment step  - not done by default

    Trimming options
      --clip_r1 [int]               Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads)
      --clip_r2 [int]               Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only)
      --three_prime_clip_r1 [int]   Instructs Trim Galore to remove bp from the 3' end of read 1 AFTER adapter/quality trimming has been performed
      --three_prime_clip_r2 [int]   Instructs Trim Galore to re move bp from the 3' end of read 2 AFTER adapter/quality trimming has been performed

    Presets:
      --pico                        Sets trimming and standedness settings for the SMARTer Stranded Total RNA-Seq Kit - Pico Input kit. Equivalent to: --forward_stranded --clip_r1 3 --three_prime_clip_r2 3
      --fcExtraAttributes           Define which extra parameters should also be included in featureCounts (default: gene_names)

    Other options:
      --outdir                      The output directory where the results will be saved
      -w/--work-dir                 The temporary directory where intermediate data will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --sampleLevel                 Used to turn of the edgeR MDS and heatmap. Set automatically when running on fewer than 3 samples
      --clusterOptions              Extra SLURM options, used in conjunction with Uppmax.config
      --maxMultiqcEmailFileSize     Theshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
      --seqCenter                   Add sequencing center in @RG line of output BAM header

    QC options:
      --skip_qc                     Skip all QC steps apart from MultiQC
      --skip_fastqc                 Skip FastQC
      --skip_rseqc                  Skip RSeQC
      --skip_genebody_coverage      Skip calculating genebody coverage
      --skip_preseq                 Skip Preseq
      --skip_dupradar               Skip dupRadar (and Picard MarkDups)
      --skip_edger                  Skip edgeR MDS plot and heatmap
      --skip_multiqc                Skip MultiQC

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}



custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

def summaryMessage() {

  log.info """=======================================================
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~\'
      |\\ | |__  __ /  ` /  \\ |__) |__         }  {
      | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                            `._,._,\'

  nf-core/rnaseq : RNA-Seq Best Practice v${workflow.manifest.version}
  ======================================================="""
  def summary = [:]
  summary['Run Name']     = custom_runName ?: workflow.runName
  summary['Reads']        = params.reads
  summary['ReadPathsFile']        = params.readPathsFile
  summary['Data Type']    = params.singleEnd ? 'Single-End' : 'Paired-End'
  summary['Genome']       = params.genome
  if( params.pico ) summary['Library Prep'] = "SMARTer Stranded Total RNA-Seq Kit - Pico Input"
  summary['Strandedness'] = ( params.unstranded ? 'None' : params.forward_stranded ? 'Forward' : params.reverse_stranded ? 'Reverse' : 'None' )
  summary['Trim R1'] = params.clip_r1
  summary['Trim R2'] = params.clip_r2
  summary["Trim 3' R1"] = params.three_prime_clip_r1
  summary["Trim 3' R2"] = params.three_prime_clip_r2
  if(params.aligner == 'star'){
      summary['Aligner'] = "STAR"
      if(params.star_index)          summary['STAR Index']   = params.star_index
      else if(params.fasta)          summary['Fasta Ref']    = params.fasta
  } else if(params.aligner == 'hisat2') {
      summary['Aligner'] = "HISAT2"
      if(params.hisat2_index)        summary['HISAT2 Index'] = params.hisat2_index
      else if(params.fasta)          summary['Fasta Ref']    = params.fasta
      if(params.splicesites)         summary['Splice Sites'] = params.splicesites
  }
  if(params.gtf)                 summary['GTF Annotation']  = params.gtf
  if(params.gff)                 summary['GFF3 Annotation']  = params.gff
  if(params.bed12)               summary['BED Annotation']  = params.bed12
  summary['Save Reference'] = params.saveReference ? 'Yes' : 'No'
  summary['Save Trimmed']   = params.saveTrimmed ? 'Yes' : 'No'
  summary['Save Intermeds'] = params.saveAlignedIntermediates ? 'Yes' : 'No'
  summary['Save Indv Quants']  = params.saveIndividualQuants ? 'Yes' : 'No'
  summary['Max Memory']     = params.max_memory
  summary['Max CPUs']       = params.max_cpus
  summary['Max Time']       = params.max_time
  summary['Output dir']     = params.outdir
  summary['Skip alignment'] = params.skip_alignment
  summary['Run salmon']     = params.run_salmon
  summary['Run exon quant'] = params.run_exon_quant
  summary['Run leafcutter'] = params.run_leafcutter
  summary['Run txrevise']   = params.run_txrevise
  summary['Working dir']    = workflow.workDir
  summary['Container']      = workflow.container
  if(workflow.revision) summary['Pipeline Release'] = workflow.revision
  summary['Current home']   = "$HOME"
  summary['Current user']   = "$USER"
  summary['Current path']   = "$PWD"
  summary['Script dir']     = workflow.projectDir
  summary['Config Profile'] = workflow.profile
  if(params.project) summary['UPPMAX Project'] = params.project
  if(params.email) {
      summary['E-mail Address'] = params.email
      summary['MultiQC maxsize'] = params.maxMultiqcEmailFileSize
  }
  log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
  log.info "========================================="


  // Show a big error message if we're running on the base config and an uppmax cluster
  if( workflow.profile == 'standard'){
      if ( "hostname".execute().text.contains('.uppmax.uu.se') ) {
          log.error "====================================================\n" +
                    "  WARNING! You are running with the default 'standard'\n" +
                    "  pipeline config profile, which runs on the head node\n" +
                    "  and assumes all software is on the PATH.\n" +
                    "  ALL JOBS ARE RUNNING LOCALLY and stuff will probably break.\n" +
                    "  Please use `-profile uppmax` to run on UPPMAX clusters.\n" +
                    "============================================================"
      }
  }
}


def emailMessage() {
    //This came from start aligment before. Right now dont have star aligner included in the workflow
    skipped_poor_alignment = []

    // Set up the e-mail variables
    def subject = "[nfcore/rnaseq] Successful: $workflow.runName"
    if(skipped_poor_alignment.size() > 0){
        subject = "[nfcore/rnaseq] Partially Successful (${skipped_poor_alignment.size()} skipped): $workflow.runName"
    }
    if(!workflow.success){
      subject = "[nfcore/rnaseq] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if(workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['skipped_poor_alignment'] = skipped_poor_alignment

    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success && !params.skip_multiqc) {
            mqc_report = multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList){
                log.warn "[nfcore/rnaseq] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
                }
        }
    } catch (all) {
        log.warn "[nfcore/rnaseq] Could not attach MultiQC report to summary email"
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.maxMultiqcEmailFileSize.toBytes() ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[nfcore/rnaseq] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nfcore/rnaseq] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Switch the embedded MIME images with base64 encoded src
    ngirnaseqlogo = new File("$baseDir/assets/nfcore-rnaseq_logo.png").bytes.encodeBase64().toString()
    email_html = email_html.replaceAll(~/cid:ngilogo/, "data:image/png;base64,$ngilogo")

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    if(skipped_poor_alignment.size() > 0){
        log.info "[nfcore/rnaseq] WARNING - ${skipped_poor_alignment.size()} samples skipped due to poor alignment scores!"
    }

    log.info "[nfcore/rnaseq] Pipeline Complete"

    if(!workflow.success){
        if( workflow.profile == 'standard'){
            if ( "hostname".execute().text.contains('.uppmax.uu.se') ) {
                log.error "====================================================\n" +
                        "  WARNING! You are running with the default 'standard'\n" +
                        "  pipeline config profile, which runs on the head node\n" +
                        "  and assumes all software is on the PATH.\n" +
                        "  This is probably why everything broke.\n" +
                        "  Please use `-profile uppmax` to run on UPPMAX clusters.\n" +
                        "============================================================"
            }
        }
    }
}