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
    file "*trimming_report.txt", emit: trimgalore_results
    file "*_fastqc.{zip,html}", emit: trimgalore_fastqc_reports
    file "where_are_my_files.txt"


    script:
    c_r1 = clip_r1 > 0 ? "--clip_r1 ${clip_r1}" : ''
    c_r2 = clip_r2 > 0 ? "--clip_r2 ${clip_r2}" : ''
    tpc_r1 = three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${three_prime_clip_r1}" : ''
    tpc_r2 = three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${three_prime_clip_r2}" : ''
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
    take: 
        raw_reads
        ch_wherearemyfiles
    main:
        trim_galore_pr(raw_reads, ch_wherearemyfiles)
    emit:
        trimmed_reads = trim_galore_pr.out.trimmed_reads
}