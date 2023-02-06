python3 write_genotypes_file.py \
    --input_vcf_list data/29632_1#149_2_PFA_GRC2_v1.0-Pf_GRC2v1.bcftools_genotyped.vcf.gz data/29632_1#55_1_PFA_GRC1_v1.0-Pf_GRC1v1.bcftools_genotyped.vcf.gz \
    --output_file_name genotype_file.tsv \
    --chromKey_file data/chromKey.txt \
    --chromosome_column_name Chrom_ID \
    --locus_column_name VarPos \
    --min_total_depth 10 \
    --het_min_allele_depth 5 \
    --het_min_allele_proportion 0.10
