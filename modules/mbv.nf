process run_mbv {
    tag "${mbv_bam.simpleName}"
    publishDir "${params.outdir}/MBV", mode: 'copy'

    input:
    file mbv_bam
    file vcf

    output:
    path "${mbv_bam.simpleName}.mbv_output.txt"

    script:
    """
    samtools index $mbv_bam
    QTLtools mbv --vcf $vcf --bam $mbv_bam --out ${mbv_bam.simpleName}.mbv_output.txt
    """
}