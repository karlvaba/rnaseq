
Channel
    .fromPath(params.mbv_vcf)
    .ifEmpty { exit 1, "VCF file is not found to perform MBV: ${params.mbv_vcf}" }
    .set { mbv_vcf_ch }


process run_mbv {
    tag "${mbv_bam.simpleName}"
    publishDir "${params.outdir}/MBV", mode: 'copy'

    input:
    file mbv_bam

    output:
    path "${mbv_bam.simpleName}.mbv_output.txt"

    script:
    """
    samtools index $mbv_bam
    QTLtools mbv --vcf ${mbv_vcf_ch.collect()} --bam $mbv_bam --out ${mbv_bam.simpleName}.mbv_output.txt
    """
}