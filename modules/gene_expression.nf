process featureCounts {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'
    
    tag "${bam_featurecounts_sorted.baseName - '.sortedByName'}"
    publishDir "${params.outdir}/featureCounts", mode: 'copy',
        saveAs: {filename ->
            if (params.saveIndividualQuants && filename.indexOf("biotype_counts") > 0) "biotype_counts/$filename"
            else if (params.saveIndividualQuants && filename.indexOf("_gene.featureCounts.txt.summary") > 0) "gene_count_summaries/$filename"
            else if (params.saveIndividualQuants && filename.indexOf("_gene.featureCounts.txt") > 0) "gene_counts/$filename"
            else null
        }

    input:
    file bam_featurecounts_sorted
    file gtf
    file biotypes_header

    output:
    path "${sample_name}_gene.featureCounts.txt", emit: geneCounts
    path "${sample_name}_gene.featureCounts.txt.summary", emit: featureCounts_logs
    path "${sample_name}_biotype_counts*mqc.{txt,tsv}", emit: featureCounts_biotype

    script:
    def featureCounts_direction = 0
    def extraAttributes = params.fcExtraAttributes ? "--extraAttributes ${params.fcExtraAttributes}" : ''
    if (forward_stranded && !unstranded) {
        featureCounts_direction = 1
    } else if (reverse_stranded && !unstranded){
        featureCounts_direction = 2
    }
    // Try to get real sample name
    sample_name = bam_featurecounts_sorted.baseName - 'ByName'
    """
    mv $bam_featurecounts_sorted ${sample_name}.bam
    featureCounts -a $gtf -g gene_id --donotsort -o ${sample_name}_gene.featureCounts.txt $extraAttributes -p -s $featureCounts_direction ${sample_name}.bam
    featureCounts -a $gtf -g gene_type --donotsort -o ${sample_name}_biotype.featureCounts.txt -p -s $featureCounts_direction ${sample_name}.bam
    cut -f 1,7 ${sample_name}_biotype.featureCounts.txt | tail -n +3 | cat $biotypes_header - >> ${sample_name}_biotype_counts_mqc.txt
    mqc_features_stat.py ${sample_name}_biotype_counts_mqc.txt -s $sample_name -f rRNA -o ${sample_name}_biotype_counts_gs_mqc.tsv
    """
}

process merge_featureCounts {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    tag "merge ${input_files.size()} files"
    publishDir "${params.outdir}/featureCounts", mode: 'copy'

    input:
    file input_files

    output:
    path 'merged_gene_counts.txt'

    script:
    //if we only have 1 file, just use cat and pipe output to csvtk. Else join all files first, and then remove unwanted column names.
    def single = input_files instanceof Path ? 1 : input_files.size()
    def merge = (single == 1) ? 'cat' : 'csvtk join -t -f "Geneid,Start,Length,End,Chr,Strand,gene_name"'
    """
    $merge $input_files | csvtk cut -t -f "-Start,-Chr,-End,-Length,-Strand" | sed 's/Aligned.sortedByCoord.out.markDups.bam//g' | sed 's/.sorted.bam//g' | csvtk rename -t -f Geneid -n phenotype_id | csvtk cut -t -f "-gene_name" > merged_gene_counts.txt
    """
}

workflow gene_expression {
    take: 
        bam_sorted 
        gtf 
        ch_biotypes_header
    main:
        featureCounts(bam_sorted, gtf, ch_biotypes_header)
        merge_featureCounts(geneCounts.toSortedList())
}