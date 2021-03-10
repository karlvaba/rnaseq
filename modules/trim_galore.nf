/*
    CHANNELS SETUP
*/


//Create a channel for input read files
if(params.readPathsFile){
    if(params.singleEnd){
        Channel.fromPath(params.readPathsFile)
        .ifEmpty { error "Cannot find any readPathsFile file in: ${params.readPathsFile}" }
        .splitCsv(header: false, sep: '\t', strip: true)
        .map{row -> [ row[0], [ file(row[1]) ] ]}
        .set { raw_reads }
    } else {
        Channel.fromPath(params.readPathsFile)
        .ifEmpty { error "Cannot find any readPathsFile file in: ${params.readPathsFile}" }
        .splitCsv(header: false, sep: '\t', strip: true)
        .map{row -> [ row[0], [ file(row[1]) , file(row[2]) ] ]}
        .set { raw_reads }
    }
} 
else if(readPaths){
    if(params.singleEnd){
        Channel
            .from(readPaths)
            .map { row -> [ row[0], [file(row[1][0])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .set { raw_reads }
    } else {
        Channel
            .from(readPaths)
            .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .set { raw_reads }
    }
} else {
    Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
        .set { raw_reads }
}

//Channel for file locations
ch_wherearemyfiles = Channel.fromPath("$baseDir/assets/where_are_my_files.txt")



/*
    TRIMGALORE PROCESSES
*/
process trim_galore_pr {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    tag "$name"
    publishDir "${params.outdir}/trim_galore", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
            else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
            else if (!params.saveTrimmed && filename == "where_are_my_files.txt") filename
            else if (params.saveTrimmed && filename != "where_are_my_files.txt") filename
            else null
        }

    input:
    tuple val(name), file(reads)
    file wherearemyfiles

    output:
    tuple val(name), file("*fq.gz"), emit: trimmed_reads
    path "*trimming_report.txt", emit: trimgalore_results
    path "*_fastqc.{zip,html}", emit: trimgalore_fastqc_reports
    path "where_are_my_files.txt"


    script:
    c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
    c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
    tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
    tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
    if (params.singleEnd) {
        """
        trim_galore --fastqc --gzip $c_r1 $tpc_r1 $reads
        """
    } else {
        """
        trim_galore --paired --fastqc --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 $reads
        """
    }
}

workflow trim_galore {
    main:
        trim_galore_pr(raw_reads, ch_wherearemyfiles)
    emit:
        trimmed_reads = trim_galore_pr.out.trimmed_reads
}