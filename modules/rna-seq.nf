/*
 * STEP 1 - FastQC
 */
process fastqc {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    when:
    !params.skip_qc && !params.skip_fastqc

    input:
    tuple val(name), file(reads)

    output:
    file "*_fastqc.{zip,html}", emit: fastqc_results

    script:
    """
    fastqc -q $reads
    """
}


/*
 * STEP 2 - Trim Galore!
 */
process trim_galore {
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
    set val(name), file("*fq.gz"), emit: trimmed_reads, trimmed_reads_salmon
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


/*
 * STEP 3 - align with STAR
 */

process star {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    tag "$prefix"
    publishDir "${params.outdir}/STAR", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".bam") == -1) "logs/$filename"
            else if (!params.saveAlignedIntermediates && filename == "where_are_my_files.txt") filename
            else if (params.saveAlignedIntermediates && filename != "where_are_my_files.txt") filename
            else null
        }

    input:
    // TODO (Nurlan Kerimov):  change the prefix to samplename in the future (did not do it because there is no test environment for changes)
    tuple samplename, file(reads)
    file index
    file gtf
    file wherearemyfiles

    output:
    // Add leafcutter and MBV bam channels to STAR alignment (here) too
    tuple file("*Log.final.out"), file ('*.bam'), emit: star_aligned
    file "*.out", emit: alignment_logs
    file "*SJ.out.tab"
    file "*Log.out", emit: star_log
    file "where_are_my_files.txt"
    file "${prefix}Aligned.sortedByCoord.out.bam.bai", emit: bam_index_rseqc, bam_index_genebody

    script:
    prefix = reads[0].toString() - ~/(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
    def star_mem = task.memory ?: params.star_memory ?: false
    def avail_mem = star_mem ? "--limitBAMsortRAM ${star_mem.toBytes() - 100000000}" : ''
    seqCenter = params.seqCenter ? "--outSAMattrRGline ID:$prefix 'CN:$params.seqCenter'" : ''
    """
    STAR --genomeDir $index \\
        --sjdbGTFfile $gtf \\
        --readFilesIn $reads  \\
        --runThreadN ${task.cpus} \\
        --twopassMode Basic \\
        --outWigType bedGraph \\
        --outSAMtype BAM SortedByCoordinate $avail_mem \\
        --readFilesCommand zcat \\
        --runDirPerm All_RWX \\
            --outFileNamePrefix $prefix $seqCenter
        
    samtools index ${prefix}Aligned.sortedByCoord.out.bam
    """
}


/*
 * STEP 3 - align with HISAT2
 */

process hisat2Align {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    tag "$samplename"
    publishDir "${params.outdir}/HISAT2", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".hisat2_summary.txt") > 0) "logs/$filename"
            else if (!params.saveAlignedIntermediates && filename == "where_are_my_files.txt") filename
            else if (params.saveAlignedIntermediates && filename != "where_are_my_files.txt") filename
            else null
        }

    input:
    set samplename, file(reads)
    file hs2_indices
    file alignment_splicesites
    file wherearemyfiles

    output:
    file "${samplename}.bam", emit: hisat2_bam
    file "${samplename}.hisat2_summary.txt", emit: alignment_logs
    file "where_are_my_files.txt"

    script:
    index_base = hs2_indices[0].toString() - ~/.\d.ht2/
    //prefix = reads[0].toString() - ~/(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
    seqCenter = params.seqCenter ? "--rg-id ${samplename} --rg CN:${params.seqCenter.replaceAll('\\s','_')}" : ''
    def rnastrandness = ''
    if (forward_stranded && !unstranded){
        rnastrandness = params.singleEnd ? '--rna-strandness F' : '--rna-strandness FR'
    } else if (reverse_stranded && !unstranded){
        rnastrandness = params.singleEnd ? '--rna-strandness R' : '--rna-strandness RF'
    }
    if (params.singleEnd) {
        """
        hisat2 -x $index_base \\
                -U $reads \\
                $rnastrandness \\
                --known-splicesite-infile $alignment_splicesites \\
                -p ${task.cpus} \\
                --met-stderr \\
                --new-summary \\
                --summary-file ${samplename}.hisat2_summary.txt $seqCenter \\
                | samtools view -bS -F 4 -F 256 - > ${samplename}.bam
        """
    } else {
        """
        hisat2 -x $index_base \\
                -1 ${reads[0]} \\
                -2 ${reads[1]} \\
                $rnastrandness \\
                --known-splicesite-infile $alignment_splicesites \\
                --no-mixed \\
                --no-discordant \\
                -p ${task.cpus} \\
                --met-stderr \\
                --new-summary \\
                --summary-file ${samplename}.hisat2_summary.txt $seqCenter \\
                | samtools view -bS -F 4 -F 8 -F 256 - > ${samplename}.bam
        """
    }
}

process hisat2_sortOutput {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    tag "${hisat2_bam.baseName}"
    publishDir "${params.outdir}/HISAT2", mode: 'copy',
        saveAs: { filename ->
            if (!params.saveAlignedIntermediates && filename == "where_are_my_files.txt") filename
            else if (params.saveAlignedIntermediates && filename != "where_are_my_files.txt") "aligned_sorted/$filename"
            else null
        }

    input:
    file hisat2_bam
    file wherearemyfiles

    output:
    file "${hisat2_bam.baseName}.sorted.bam", emit: bam_count, bam_rseqc, bam_preseq, bam_markduplicates, bam_featurecounts, bam_for_genebody, leafcutter_bam, mbv_bam
    file "${hisat2_bam.baseName}.sorted.bam.bai", emit bam_index_rseqc, bam_index_genebody
    file "where_are_my_files.txt"

    script:
    def avail_mem = task.memory ? "-m ${task.memory.toBytes() / task.cpus}" : ''
    """
    samtools sort \\
        $hisat2_bam \\
        -@ ${task.cpus} $avail_mem \\
        -o ${hisat2_bam.baseName}.sorted.bam
    samtools index ${hisat2_bam.baseName}.sorted.bam
    """
}


process sort_by_name_BAM {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    tag "${bam_featurecounts.baseName - '.sorted'}"

    input:
    file bam_featurecounts

    output:
    file "${bam_featurecounts.baseName}ByName.bam", emit: bam_featurecounts_sorted, bam_count_exons

    script:
    def avail_mem = task.memory ? "-m ${task.memory.toBytes() / (task.cpus + 2)}" : ''
    """
    samtools sort -n \\
        $bam_featurecounts \\
        -@ ${task.cpus} $avail_mem \\
        -o ${bam_featurecounts.baseName}ByName.bam
    """
}



/*
 * STEP 3Salmon.1 - quant transcripts with Salmon
 */

process salmon_quant {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    tag "$samplename - ${index.baseName}"
    publishDir "${params.outdir}/Salmon/quant/${index.baseName}/", mode: 'copy',
        saveAs: {filename -> if (params.saveIndividualQuants && filename.indexOf(".quant.sf") > 0) filename else null }

    input:
    tuple samplename, file(reads)
    each index

    output:
    tuple val(index.baseName), file("${samplename}.quant.edited.sf"), emit: salmon_merge_tx_ch
    file '*.quant.sf'
    
    script:
    def strandedness = params.unstranded ? 'U' : 'SR'
    if (params.singleEnd) {
        """
        salmon quant --seqBias --useVBOpt --gcBias \\
                        --libType $strandedness \\
                        --index ${index} \\
                        -r ${reads[0]} \\
                        -p ${task.cpus} \\
                        -o .
        mv quant.sf ${samplename}.quant.sf
        cat ${samplename}.quant.sf | csvtk cut -t -f "-Length,-EffectiveLength" | sed '1s/TPM/${samplename}_TPM/g' | sed '1s/NumReads/${samplename}_NumReads/g' > ${samplename}.quant.edited.sf
        """
    } else {
        """
        salmon quant --seqBias --useVBOpt --gcBias \\
                        --libType I$strandedness \\
                        --index $index \\
                        -1 ${reads[0]} \\
                        -2 ${reads[1]} \\
                        -p ${task.cpus} \\
                        -o .
        mv quant.sf ${samplename}.quant.sf
        cat ${samplename}.quant.sf | csvtk cut -t -f "-Length,-EffectiveLength" | sed '1s/TPM/${samplename}_TPM/g' | sed '1s/NumReads/${samplename}_NumReads/g' > ${samplename}.quant.edited.sf
        """
    }
}


/*
 * STEP 3Salmon.2 - merge salmon outputs
 */

process salmon_merge {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    tag "merge_salmon_${index}"
    publishDir "${params.outdir}/Salmon/merged_counts/", mode: 'copy',
        saveAs: {filename -> if (filename.indexOf("TPM.merged.txt") > 0) "TPM/$filename"
        else if (filename.indexOf(".NumReads.merged.txt") > 0) "NumReads/$filename"
        else null }

    input:
    tuple index, file(input_files)

    output:
    file '*merged.txt'

    script:
    //if we only have 1 file, just use cat and pipe output to csvtk. Else join all files first, and then remove unwanted column names.
    def single = input_files instanceof Path ? 1 : input_files.size()
    def merge = (single == 1) ? 'cat' : 'csvtk join -t -f "Name"'
    """
    $merge $input_files | csvtk rename -t -f Name -n phenotype_id > merged_TPMS_NumReads.tsv
    csvtk cut -t -F -f -"*_NumReads" merged_TPMS_NumReads.tsv | sed 's/_TPM//g' > ${index}.TPM.merged.txt
    csvtk cut -t -F -f -"*_TPM" merged_TPMS_NumReads.tsv | sed 's/_NumReads//g' > ${index}.NumReads.merged.txt
    """
}


/*
 * Leafcutter quantification preparation step
 */

process leafcutter_bam_to_junc {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    tag "${leafcutter_bam.baseName}"
    
    input:
    file leafcutter_bam

    output:
    file '*.junc', emit: leafcutter_junc
    
    script:
    """
    samtools view $leafcutter_bam | python $baseDir/bin/leafcutter/filter_cs.py | $baseDir/bin/leafcutter/sam2bed.pl --use-RNA-strand - ${leafcutter_bam.baseName}.bed
    $baseDir/bin/leafcutter/bed2junc.pl ${leafcutter_bam.baseName}.bed ${leafcutter_bam.baseName}.junc
    """
}


/*
 * Leafcutter quantification step
 */

process leafcutter_cluster_junctions {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    tag "${junc_files.baseName}"
    publishDir "${params.outdir}/leafcutter", mode: 'copy'

    input:
    file junc_files

    output:
    file 'leafcutter_perind*.formatted.gz'
    file 'junction_files.tar.gz'
    
    script:
    """
    tar -czf junction_files.tar.gz -T $junc_files
    python $baseDir/bin/leafcutter/leafcutter_cluster.py -j $junc_files -r . -m ${params.leafcutter_min_split_reads} -l ${params.leafcutter_max_intron_length}
    zcat leafcutter_perind.counts.gz | csvtk space2tab | csvtk rename -t -f chrom -n phenotype_id | sed 's/.sorted//g' | gzip -c > leafcutter_perind.counts.formatted.gz 
    zcat leafcutter_perind_numers.counts.gz | sed '1s/^/phenotype_id /' | sed 's/.sorted//g' | csvtk space2tab | gzip -c > leafcutter_perind_numers.counts.formatted.gz
    """        
}


/*
 * Quantify exon expression - featureCounts (exon level)
 */

process count_exons {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    tag "${bam.simpleName}"
    
    publishDir "${params.outdir}/dexseq_exon_counts/quant_files", mode: 'copy',
        saveAs: {filename -> if (params.saveIndividualQuants && filename.indexOf(".exoncount.txt") > 0) filename else null }


    input:
    file bam
    file gff

    output:
    file "${bam.simpleName}.exoncount.txt", emit: merge_exon_count_ch

    script:
    def featureCounts_direction = 0
    if (forward_stranded && !unstranded) {
        featureCounts_direction = 1
    } else if (reverse_stranded && !unstranded){
        featureCounts_direction = 2
    }
    """
    featureCounts -p -t exonic_part -s $featureCounts_direction -f -O -a $gff -o ${bam.simpleName}.exoncount.txt $bam
    """
}


/*
 * Merge exon counts files
 */

process exon_count_merge {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    tag "merge exon ${input_files.size()} files"
    publishDir "${params.outdir}/dexseq_exon_counts", mode: 'copy'

    input:
    file input_files

    output:
    file 'merged_exon_counts.tsv'

    script:
    //if we only have 1 file, just use cat and pipe output to csvtk. Else join all files first, and then remove unwanted column names.
    def single = input_files instanceof Path ? 1 : input_files.size()
    def merge = (single == 1) ? 'cat' : 'csvtk join -t -f "Geneid,Start,Length,End,Chr,Strand"'
    """
    $merge $input_files | sed 's/.sortedByName.bam//g' | awk '\$1=\$1"_"\$2"_"\$3"_"\$4' OFS='\t' | csvtk rename -t -f Geneid_Chr_Start_End -n phenotype_id | csvtk cut -t -f "-Chr,-Start,-End,-Strand,-Length" > merged_exon_counts.tsv
    """
}


/*
 * Run QTLtools MBV
 */

process run_mbv {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    tag "${mbv_bam.simpleName}"
    publishDir "${params.outdir}/MBV", mode: 'copy'

    input:
    file mbv_bam
    file vcf

    output:
    file "${mbv_bam.simpleName}.mbv_output.txt"

    script:
    """
    samtools index $mbv_bam
    QTLtools mbv --vcf $vcf --bam $mbv_bam --out ${mbv_bam.simpleName}.mbv_output.txt
    """
}


/*
 * Parse software version numbers
 */
process get_software_versions {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    output:
    file 'software_versions_mqc.yaml', emit: software_versions_yaml

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
    STEP 4 - RSeQC analysis
*/
process rseqc {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'
    
    tag "${bam_rseqc.baseName - '.sorted'}"
    publishDir "${params.outdir}/rseqc" , mode: 'copy',
        saveAs: {filename ->
                    if (filename.indexOf("bam_stat.txt") > 0)                      "bam_stat/$filename"
            else if (filename.indexOf("infer_experiment.txt") > 0)              "infer_experiment/$filename"
            else if (filename.indexOf("read_distribution.txt") > 0)             "read_distribution/$filename"
            else if (filename.indexOf("read_duplication.DupRate_plot.pdf") > 0) "read_duplication/$filename"
            else if (filename.indexOf("read_duplication.DupRate_plot.r") > 0)   "read_duplication/rscripts/$filename"
            else if (filename.indexOf("read_duplication.pos.DupRate.xls") > 0)  "read_duplication/dup_pos/$filename"
            else if (filename.indexOf("read_duplication.seq.DupRate.xls") > 0)  "read_duplication/dup_seq/$filename"
            else if (filename.indexOf("RPKM_saturation.eRPKM.xls") > 0)         "RPKM_saturation/rpkm/$filename"
            else if (filename.indexOf("RPKM_saturation.rawCount.xls") > 0)      "RPKM_saturation/counts/$filename"
            else if (filename.indexOf("RPKM_saturation.saturation.pdf") > 0)    "RPKM_saturation/$filename"
            else if (filename.indexOf("RPKM_saturation.saturation.r") > 0)      "RPKM_saturation/rscripts/$filename"
            else if (filename.indexOf("inner_distance.txt") > 0)                "inner_distance/$filename"
            else if (filename.indexOf("inner_distance_freq.txt") > 0)           "inner_distance/data/$filename"
            else if (filename.indexOf("inner_distance_plot.r") > 0)             "inner_distance/rscripts/$filename"
            else if (filename.indexOf("inner_distance_plot.pdf") > 0)           "inner_distance/plots/$filename"
            else if (filename.indexOf("junction_plot.r") > 0)                   "junction_annotation/rscripts/$filename"
            else if (filename.indexOf("junction.xls") > 0)                      "junction_annotation/data/$filename"
            else if (filename.indexOf("splice_events.pdf") > 0)                 "junction_annotation/events/$filename"
            else if (filename.indexOf("splice_junction.pdf") > 0)               "junction_annotation/junctions/$filename"
            else if (filename.indexOf("junctionSaturation_plot.pdf") > 0)       "junction_saturation/$filename"
            else if (filename.indexOf("junctionSaturation_plot.r") > 0)         "junction_saturation/rscripts/$filename"
            else filename
        }

    when:
    !params.skip_qc && !params.skip_rseqc

    input:
    file bam_rseqc
    file index
    file bed12

    output:
    file "*.{txt,pdf,r,xls}", emit: rseqc_results

    script:
    """
    infer_experiment.py -i $bam_rseqc -r $bed12 > ${bam_rseqc.baseName}.infer_experiment.txt
    junction_annotation.py -i $bam_rseqc -o ${bam_rseqc.baseName}.rseqc -r $bed12
    bam_stat.py -i $bam_rseqc 2> ${bam_rseqc.baseName}.bam_stat.txt
    junction_saturation.py -i $bam_rseqc -o ${bam_rseqc.baseName}.rseqc -r $bed12 2> ${bam_rseqc.baseName}.junction_annotation_log.txt
    inner_distance.py -i $bam_rseqc -o ${bam_rseqc.baseName}.rseqc -r $bed12
    read_distribution.py -i $bam_rseqc -r $bed12 > ${bam_rseqc.baseName}.read_distribution.txt
    read_duplication.py -i $bam_rseqc -o ${bam_rseqc.baseName}.read_duplication
    """
}


/*
* Step 4.1 Rseqc create BigWig coverage
*/

process createBigWig {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'
    
    tag "${bam.baseName - 'sortedByCoord.out'}"
    publishDir "${params.outdir}/bigwig", mode: 'copy'

    when:
    !params.skip_qc && !params.skip_genebody_coverage

    input:
    file bam
    file index

    output:
    file "*.bigwig", emit: bigwig_for_genebody

    script:
    """
    bamCoverage -b $bam -p ${task.cpus} -o ${bam.baseName}.bigwig
    """
}
/*
    * Step 4.2 Rseqc genebody_coverage
    */
process genebody_coverage {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    tag "${bigwig.baseName}"
        publishDir "${params.outdir}/rseqc" , mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("geneBodyCoverage.curves.pdf") > 0)       "geneBodyCoverage/$filename"
            else if (filename.indexOf("geneBodyCoverage.r") > 0)           "geneBodyCoverage/rscripts/$filename"
            else if (filename.indexOf("geneBodyCoverage.txt") > 0)         "geneBodyCoverage/data/$filename"
            else if (filename.indexOf("log.txt") > -1) false
            else filename
        }

    when:
    !params.skip_qc && !params.skip_genebody_coverage

    input:
    file bigwig
    file bed12

    output:
    file "*.{txt,pdf,r}", emit: genebody_coverage_results

    script:
    """
    geneBody_coverage2.py \\
        -i $bigwig \\
        -o ${bigwig.baseName}.rseqc.txt \\
        -r $bed12
    """
}

/*
    * STEP 5 - preseq analysis
    */
process preseq {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    tag "${bam_preseq.baseName - '.sorted'}"
    publishDir "${params.outdir}/preseq", mode: 'copy'

    when:
    !params.skip_qc && !params.skip_preseq

    input:
    file bam_preseq

    output:
    file "${bam_preseq.baseName}.ccurve.txt", emit: preseq_results

    script:
    """
    preseq lc_extrap -v -B $bam_preseq -o ${bam_preseq.baseName}.ccurve.txt
    """
}


/*
    * STEP 6 Mark duplicates
    */
process markDuplicates {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    tag "${bam.baseName - '.sorted'}"
    publishDir "${params.outdir}/markDuplicates", mode: 'copy',
        saveAs: {filename -> filename.indexOf("_metrics.txt") > 0 ? "metrics/$filename" : "$filename"}

    when:
    !params.skip_qc && !params.skip_dupradar

    input:
    file bam

    output:
    file "${bam.baseName}.markDups.bam", emit: bam_md
    file "${bam.baseName}.markDups_metrics.txt", emit: picard_results
    file "${bam.baseName}.markDups.bam.bai"

    script:
    if( !task.memory ){
        log.info "[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this."
        avail_mem = 3
    } else {
        avail_mem = task.memory.toGiga()
    }
    """
    picard -Xmx${avail_mem}g MarkDuplicates \\
        INPUT=$bam \\
        OUTPUT=${bam.baseName}.markDups.bam \\
        METRICS_FILE=${bam.baseName}.markDups_metrics.txt \\
        REMOVE_DUPLICATES=false \\
        ASSUME_SORTED=true \\
        PROGRAM_RECORD_ID='null' \\
        VALIDATION_STRINGENCY=LENIENT
    samtools index ${bam.baseName}.markDups.bam
    """
}


/*
    * STEP 7 - dupRadar
    */
process dupradar {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    tag "${bam_md.baseName - '.sorted.markDups'}"
    publishDir "${params.outdir}/dupradar", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_duprateExpDens.pdf") > 0) "scatter_plots/$filename"
            else if (filename.indexOf("_duprateExpBoxplot.pdf") > 0) "box_plots/$filename"
            else if (filename.indexOf("_expressionHist.pdf") > 0) "histograms/$filename"
            else if (filename.indexOf("_dupMatrix.txt") > 0) "gene_data/$filename"
            else if (filename.indexOf("_duprateExpDensCurve.txt") > 0) "scatter_curve_data/$filename"
            else if (filename.indexOf("_intercept_slope.txt") > 0) "intercepts_slopes/$filename"
            else "$filename"
        }

    when:
    !params.skip_qc && !params.skip_dupradar

    input:
    file bam_md
    file gtf

    output:
    file "*.{pdf,txt}", emit: dupradar_results

    script: // This script is bundled with the pipeline, in nfcore/rnaseq/bin/
    def dupradar_direction = 0
    if (forward_stranded && !unstranded) {
        dupradar_direction = 1
    } else if (reverse_stranded && !unstranded){
        dupradar_direction = 2
    }
    def paired = params.singleEnd ? 'single' :  'paired'
    """
    dupRadar.r $bam_md $gtf $dupradar_direction $paired ${task.cpus}
    """
}

/*
    * STEP 8 Feature counts
    */
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
    file "${sample_name}_gene.featureCounts.txt", emit: geneCounts, featureCounts_to_merge
    file "${sample_name}_gene.featureCounts.txt.summary", emit: featureCounts_logs
    file "${sample_name}_biotype_counts*mqc.{txt,tsv}", emit: featureCounts_biotype

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

/*
    * STEP 9 - Merge featurecounts
    */
process merge_featureCounts {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    tag "merge ${input_files.size()} files"
    publishDir "${params.outdir}/featureCounts", mode: 'copy'

    input:
    file input_files

    output:
    file 'merged_gene_counts.txt'

    script:
    //if we only have 1 file, just use cat and pipe output to csvtk. Else join all files first, and then remove unwanted column names.
    def single = input_files instanceof Path ? 1 : input_files.size()
    def merge = (single == 1) ? 'cat' : 'csvtk join -t -f "Geneid,Start,Length,End,Chr,Strand,gene_name"'
    """
    $merge $input_files | csvtk cut -t -f "-Start,-Chr,-End,-Length,-Strand" | sed 's/Aligned.sortedByCoord.out.markDups.bam//g' | sed 's/.sorted.bam//g' | csvtk rename -t -f Geneid -n phenotype_id | csvtk cut -t -f "-gene_name" > merged_gene_counts.txt
    """
}

/*
    * STEP 11 - edgeR MDS and heatmap
    */
process sample_correlation {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    tag "${input_files[0].toString() - '.sorted_gene.featureCounts.txt' - 'Aligned'}"
    publishDir "${params.outdir}/sample_correlation", mode: 'copy'

    when:
    !params.skip_qc && !params.skip_edger

    input:
    file input_files
    val num_bams
    file mdsplot_header
    file heatmap_header

    output:
    file "*.{txt,pdf,csv}", emit: sample_correlation_results

    when:
    num_bams > 2 && (!params.sampleLevel)

    script: // This script is bundled with the pipeline, in nfcore/rnaseq/bin/
    """
    edgeR_heatmap_MDS.r $input_files
    cat $mdsplot_header edgeR_MDS_Aplot_coordinates_mqc.csv >> tmp_file
    mv tmp_file edgeR_MDS_Aplot_coordinates_mqc.csv
    cat $heatmap_header log2CPM_sample_distances_mqc.csv >> tmp_file
    mv tmp_file log2CPM_sample_distances_mqc.csv
    """
}

/*
    * Pipeline parameters to go into MultiQC report
    */
process workflow_summary_mqc {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    when:
    !params.skip_multiqc

    output:
    file 'workflow_summary_mqc.yaml', emit: workflow_summary_yaml

    exec:
    def yaml_file = task.workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nfcore-rnaseq-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nfcore/rnaseq Workflow Summary'
    section_href: 'https://github.com/nf-core/rnaseq'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()
}

/*
    * STEP 12 MultiQC
    */
process multiqc {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    when:
    !params.skip_multiqc

    input:
    file multiqc_config
    file (fastqc:'fastqc/*')
    file ('trimgalore/*')
    file ('alignment/*')
    file ('rseqc/*')
    file ('rseqc/*')
    file ('preseq/*')
    file ('dupradar/*')
    file ('featureCounts/*')
    file ('featureCounts_biotype/*')
    file ('sample_correlation_results/*') // If the Edge-R is not run create an Empty array
    file ('software_versions/*')
    file ('workflow_summary/*')

    output:
    file "*multiqc_report.html", emit: multiqc_report
    file "*_data"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc . -f $rtitle $rfilename --config $multiqc_config \\
        -m custom_content -m picard -m preseq -m rseqc -m featureCounts -m hisat2 -m star -m cutadapt -m fastqc
    """
}


/*
 * STEP 13 - Output Description HTML
 */
process output_documentation {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    file output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}