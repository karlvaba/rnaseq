nextflow run main.nf\
 -profile tartu_hpc,eqtl_catalogue\
 --readPathsFile /gpfs/space/projects/GTEx/GTExV8_group_readpaths/group_8_readpaths.tsv\
 --unstranded\
 --aligner 'hisat2'\
 --run_tx_exp_quant\
 --run_txrevise\
 --run_exon_quant\
 --saveReference\
 --run_mbv\
 --mbv_vcf /gpfs/space/projects/GTEx/genotypes/processed/GTEx.MAF001.vcf.gz\
 --outdir /gpfs/space/projects/GTEx/rnaseq/group_8\
 -resume
