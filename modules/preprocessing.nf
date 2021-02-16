/*
 * PREPROCESSING - Build STAR index
 */

process makeSTARindex {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    tag "$fasta"
    publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                saveAs: { params.saveReference ? it : null }, mode: 'copy'

    input:
    file fasta
    file gtf

    output:
    file "star", emit: star_index

    script:
    def avail_mem = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    """
    mkdir star
    STAR \\
        --runMode genomeGenerate \\
        --runThreadN ${task.cpus} \\
        --sjdbGTFfile $gtf \\
        --genomeDir star/ \\
        --genomeFastaFiles $fasta \\
        $avail_mem
    """
}

/*
 * PREPROCESSING - Build HISAT2 splice sites file
 */

process makeHisatSplicesites {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    tag "$gtf"
    publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                saveAs: { params.saveReference ? it : null }, mode: 'copy'

    input:
    file gtf

    output:
    file "${gtf.baseName}.hisat2_splice_sites.txt", emit: indexing_splicesites, alignment_splicesites //TODO: does it work like this?

    script:
    """
    hisat2_extract_splice_sites.py $gtf > ${gtf.baseName}.hisat2_splice_sites.txt
    """
}

/*
 * PREPROCESSING - Build HISAT2 index
 */
process makeHISATindex {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1' 

    tag "$fasta"
    publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                saveAs: { params.saveReference ? it : null }, mode: 'copy'

    input:
    file fasta
    file indexing_splicesites
    file gtf

    output:
    file "${fasta.baseName}.*.ht2", emit: hs2_indices

    script:
    if( !task.memory ){
        log.info "[HISAT2 index build] Available memory not known - defaulting to 0. Specify process memory requirements to change this."
        avail_mem = 0
    } else {
        log.info "[HISAT2 index build] Available memory: ${task.memory}"
        avail_mem = task.memory.toGiga()
    }
    if( avail_mem > params.hisatBuildMemory ){
        log.info "[HISAT2 index build] Over ${params.hisatBuildMemory} GB available, so using splice sites and exons in HISAT2 index"
        extract_exons = "hisat2_extract_exons.py $gtf > ${gtf.baseName}.hisat2_exons.txt"
        ss = "--ss $indexing_splicesites"
        exon = "--exon ${gtf.baseName}.hisat2_exons.txt"
    } else {
        log.info "[HISAT2 index build] Less than ${params.hisatBuildMemory} GB available, so NOT using splice sites and exons in HISAT2 index."
        log.info "[HISAT2 index build] Use --hisatBuildMemory [small number] to skip this check."
        extract_exons = ''
        ss = ''
        exon = ''
    }
    """
    $extract_exons
    hisat2-build -p ${task.cpus} $ss $exon $fasta ${fasta.baseName}.hisat2_index
    """
}


/*
 * PREPROCESSING - txrevise gff3 to fasta
 */
process gff_to_fasta {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    tag "${txrevise_gff.baseName}"
    publishDir path: { params.saveReference ? "${params.outdir}/Salmon/salmon_fasta" : params.outdir },
                saveAs: { params.saveReference ? it : null }, mode: 'copy'

    input:
    file txrevise_gff
    file genome_fasta

    output:
    file "${txrevise_gff.baseName}.fa", emit: txrevise_fasta_ch
    
    script:
    """
    gffread -w ${txrevise_gff.baseName}.fa -g $genome_fasta $txrevise_gff
    """
}

/*
 * PREPROCESSING - Build Salmon index
 */
process makeSalmonIndex {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    tag "${fasta.baseName}"
    publishDir path: { params.saveReference ? "${params.outdir}/Salmon/salmon_index" : params.outdir },
                saveAs: { params.saveReference ? it : null }, mode: 'copy'
    
    input:
    file fasta

    output:
    file "${fasta.baseName}.index", emit: salmon_index
    
    script:
    """
    salmon index -t ${fasta} -i ${fasta.baseName}.index
    """
}


/*
 * PREPROCESSING - Convert GFF3 to GTF
 */

process convertGFFtoGTF {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    tag "$gff"

    input:
    file gff from gffFile

    output:
    file "${gff.baseName}.gtf", emit: gtf_makeSTARindex, gtf_makeHisatSplicesites, gtf_makeHISATindex, gtf_makeBED12,
        gtf_star, gtf_dupradar, gtf_featureCounts, gtf_dexseq

    script:
    """
    gffread  $gff -T -o ${gff.baseName}.gtf
    """
}


/*
 * PREPROCESSING - Build Exon GFF for dexseq
 */

process makeDexSeqExonGFF {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    tag "${gtf.baseName}"
    publishDir path: { params.saveReference ? "${params.outdir}/dexseq_exon_counts" : params.outdir },
                saveAs: { params.saveReference ? it : null }, mode: 'copy'

    input:
    file gtf

    output:
    file "${gtf.baseName}.patched_contigs.DEXSeq.gff", emit: dexseq_gff_count_exons
    
    script:
    """
    cat $gtf | sed 's/chrM/chrMT/;s/chr//' > ${gtf.baseName}.patched_contigs.gtf
    $baseDir/bin/dexseq/dexseq_prepare_annotation.py ${gtf.baseName}.patched_contigs.gtf ${gtf.baseName}.patched_contigs.DEXSeq.gff
    """
}


/*
 * PREPROCESSING - Build BED12 file
 */

process makeBED12 {
    container = 'quay.io/eqtlcatalogue/rnaseq:v20.11.1'

    tag "$gtf"
    publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                saveAs: { params.saveReference ? it : null }, mode: 'copy'

    input:
    file gtf

    output:
    file "${gtf.baseName}.bed", emit: bed_rseqc, bed_genebody_coverage

    script: // This script is bundled with the pipeline, in nfcore/rnaseq/bin/
    """
    gtf2bed $gtf > ${gtf.baseName}.bed
    """
}
