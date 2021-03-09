//Required for txrevise. Maybe should be in some other file
process gff_to_fasta {
    tag "${txrevise_gff.baseName}"
    publishDir path: { params.saveReference ? "${params.outdir}/Salmon/salmon_fasta" : params.outdir },
                saveAs: { params.saveReference ? it : null }, mode: 'copy'

    input:
    file txrevise_gff
    file genome_fasta

    output:
    path "${txrevise_gff.baseName}.fa", emit: txrevise_fasta_ch
    
    script:
    """
    gffread -w ${txrevise_gff.baseName}.fa -g $genome_fasta $txrevise_gff
    """
}


//Salmon processes
process makeSalmonIndex {
    tag "${fasta.baseName}"
    publishDir path: { params.saveReference ? "${params.outdir}/Salmon/salmon_index" : params.outdir },
                saveAs: { params.saveReference ? it : null }, mode: 'copy'
    
    input:
    file fasta

    output:
    path "${fasta.baseName}.index", emit: salmon_index
    
    script:
    """
    salmon index -t ${fasta} -i ${fasta.baseName}.index
    """
}

process salmon_quant {
    tag "$samplename - ${index.baseName}"
    publishDir "${params.outdir}/Salmon/quant/${index.baseName}/", mode: 'copy',
        saveAs: {filename -> if (params.saveIndividualQuants && filename.indexOf(".quant.sf") > 0) filename else null }

    input:
    tuple samplename, file(reads)
    each index

    output:
    tuple val(index.baseName), file("${samplename}.quant.edited.sf"), emit: salmon_merge_tx_ch
    path '*.quant.sf'
    
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
process salmon_merge {
    tag "merge_salmon_${index}"
    publishDir "${params.outdir}/Salmon/merged_counts/", mode: 'copy',
        saveAs: {filename -> if (filename.indexOf("TPM.merged.txt") > 0) "TPM/$filename"
        else if (filename.indexOf(".NumReads.merged.txt") > 0) "NumReads/$filename"
        else null }

    input:
    tuple index, file(input_files)

    output:
    path '*merged.txt'

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


workflow transcript_expression {
    take:
        salmon_fasta_ch
        trimmed_reads
    main:
        makeSalmonIndex(salmon_fasta_ch)
        salmon_quant(trimmed_reads, makeSalmonIndex.out.salmon_index)
        salmon_merge(salmon_quant.out.salmon_merge_tx_ch.groupTuple())
}