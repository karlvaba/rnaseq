#!/usr/bin/env nextflow
/*
===============================================================
 nf-core/rnaseq
===============================================================
 RNA-Seq Analysis Pipeline. Started March 2016.
 #### Homepage / Documentation
 https://github.com/nf-core/rnaseq
 #### Authors
 Phil Ewels @ewels <phil.ewels@scilifelab.se>
 Rickard Hammarén @Hammarn  <rickard.hammaren@scilifelab.se>
---------------------------------------------------------------
*/
nextflow.enable.dsl=2

include { helpMessage; emailMessage } from './modules/utils'

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}



// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// Reference index path configuration
// Define these here - after the profiles are loaded with the iGenomes paths
params.star_index = params.genome ? params.genomes[ params.genome ].star ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.gff = params.genome ? params.genomes[ params.genome ].gff ?: false : false
params.hisat2_index = params.genome ? params.genomes[ params.genome ].hisat2 ?: false : false



// Profile validations
if( workflow.profile == 'uppmax' || workflow.profile == 'uppmax-devel' ){
    if ( !params.project ) exit 1, "No UPPMAX project ID found! Use --project"
}

//AWSBatch sanity checking
if(workflow.profile == 'awsbatch'){
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
}



/*
    Main workflow
*/

include { trim_galore } from './modules/trim_galore'
include { align_hisat2 } from './modules/hisat2'
include { transcript_expression } from './modules/transcript'
include { output_docs; get_software_versions; summaryMessage} from './modules/utils'

workflow {
    summaryMessage()
    get_software_versions()

    ch_wherearemyfiles = Channel.fromPath("$baseDir/assets/where_are_my_files.txt")
    trim_galore(ch_wherearemyfiles)

    if(!params.skip_alignment) {
        //Also runs expressions if required by params
        align_hisat2(trim_galore.out.trimmed_reads, ch_wherearemyfiles)
    }

    if(params.run_salmon || params.run_txrevise ) {
        //Also includes tx_revise if required by params
        transcript_expression(trim_galore.out.trimmed_reads)
    }

    output_docs()
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {
    emailMessage()
}
