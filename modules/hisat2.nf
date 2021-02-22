process makeHisatSplicesites {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    tag "$gtf"
    publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                saveAs: { params.saveReference ? it : null }, mode: 'copy'

    input:
    file gtf
    output:
    path "${gtf.baseName}.hisat2_splice_sites.txt", emit: splicesites

    script:
    """
    hisat2_extract_splice_sites.py $gtf > ${gtf.baseName}.hisat2_splice_sites.txt
    """
}


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
    tuple val(clip_r1), val(clip_r2), val(three_prime_clip_r1), val(three_prime_clip_r2), val(forward_stranded), val(reverse_stranded), val(unstranded)
    tuple val(samplename), file(reads)
    file hs2_indices
    file alignment_splicesites
    file wherearemyfiles

    output:
    path "${samplename}.bam", emit: hisat2_bam
    path "${samplename}.hisat2_summary.txt", emit: alignment_logs
    path "where_are_my_files.txt"

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
    path "${hisat2_bam.baseName}.sorted.bam", emit: bam
    path "${hisat2_bam.baseName}.sorted.bam.bai", emit: bam_index
    path "where_are_my_files.txt"

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
    path "${bam_featurecounts.baseName}ByName.bam", emit: bam_sorted

    script:
    def avail_mem = task.memory ? "-m ${task.memory.toBytes() / (task.cpus + 2)}" : ''
    """
    samtools sort -n \\
        $bam_featurecounts \\
        -@ ${task.cpus} $avail_mem \\
        -o ${bam_featurecounts.baseName}ByName.bam
    """
}


workflow align_hisat2 {
    take:
        variables 
        gtf
        trimmed_reads
        hs2_indices
        ch_wherearemyfiles
    main:
        makeHisatSplicesites(gtf)
        hisat2Align(variables, trimmed_reads, hs2_indices.collect(), makeHisatSplicesites.out.splicesites, ch_wherearemyfiles.collect())
        hisat2_sortOutput(hisat2Align.out.hisat2_bam, ch_wherearemyfiles.collect())
        sort_by_name_BAM(hisat2_sortOutput.out.bam)
    emit:
        bam_sorted = sort_by_name_BAM.out.bam_sorted
}