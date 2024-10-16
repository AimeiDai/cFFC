#!/bin/bash

# cd /data/users/aimei/cFF/miR_expression_embryo/fastq_file
# sra=(SRR014367.2 SRR068998.2 SRR069234.2 SRR069262.2 SRR069503.2 SRR069504.2 SRR069838.2 SRR069839.2 SRR069840.2)
# ### step 1, convert sra to fastq 
# for i in ${!sra[@]}
# #for i in $(seq 0 5)
# do
# 	mv ${sra[$i]} ${sra[$i]}.sra
# 	fastq-dump ${sra[$i]}.sra --gzip 
# done
#rm *.sra

cd /public1/users/aimei/cFF_revised/2.miR_TF_expression/miR_expression

# sra=(SRR013600 SRR013602 SRR014367 SRR015372 SRR028728 SRR028729 SRR14149905 SRR14149914 SRR14149915 SRR14149916 SRR14149917 SRR14149918)
# ### step 2, trim adapter sequences
# mkdir trimmed_fastq
# for i in ${!sra[@]}
# #for i in $(seq 0 5)
# do
# 	{
# 		cutadapt -a CTGTAGGCACCATCAATC -O 3 -m 18 -o trimmed_fastq/${sra[$i]}.fastq.gz fastq_file/${sra[$i]}_1.fastq.gz
# 	}&
# 	if [[ `expr \( $i + 1 \) % 3` -eq 0 ]]; then
# 		wait
# 	fi
# done

# sra=(SRR14507258 SRR14507259 SRR14507260 SRR14507261 SRR14507354 SRR14507355 SRR14507356 SRR14507357 SRR14507358 SRR14507359 SRR14507360 SRR14507361 SRR14507362 SRR14507364 SRR14507365 SRR14507366 SRR14507367 SRR14507368 SRR14507369 SRR14507370 SRR14507371 SRR14507372 SRR14507373)
# for i in ${!sra[@]}
# #for i in $(seq 0 5)
# do
# 	{
# 		cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAT -O 3 -m 18 -o trimmed_fastq/${sra[$i]}.fastq.gz fastq_file/${sra[$i]}_1.fastq.gz
# 	}&
# 	if [[ `expr \( $i + 1 \) % 3` -eq 0 ]]; then
# 		wait
# 	fi
# done



fq=(SRR013600 SRR013602 SRR014367 SRR015372 SRR028728 SRR028729 SRR14149905 SRR14149914 SRR14149915 SRR14149916 SRR14149917 SRR14149918 SRR14507258 SRR14507259 SRR14507260 SRR14507261 SRR14507354 SRR14507355 SRR14507356 SRR14507357 SRR14507358 SRR14507359 SRR14507360 SRR14507361 SRR14507362 SRR14507364 SRR14507365 SRR14507366 SRR14507367 SRR14507368 SRR14507369 SRR14507370 SRR14507371 SRR14507372 SRR14507373)
mkdir fastq_unaligned_tRNA
mkdir fastq_aligned_tRNA

# for i in ${!fq[@]}
# do
# 	{
# 		bowtie -x /public1/users/aimei/genomes/Drosophila_melanogaster/small_RNA/tRNA_sequence/tRNA_dme -p 4 trimmed_fastq/${fq[$i]}.fastq.gz --un fastq_unaligned_tRNA/${fq[$i]}.fastq.gz -S fastq_aligned_tRNA/${fq[$i]}.sam 2>fastq_unaligned_tRNA/${fq[$i]}.map.info
# 	}&
# 	if [[ `expr \( $i + 1 \) % 3` -eq 0 ]]; then
# 		wait
# 	fi
# done

# mkdir fastq_unaligned_TE_consensus
# mkdir fastq_aligned_TE_consensus


# for i in ${!fq[@]}
# do
# 	{
# 		bowtie -x /public1/users/aimei/genomes/Drosophila_melanogaster/small_RNA/TE_annotation/bowtie_index/consensus -p 4 fastq_unaligned_tRNA/${fq[$i]}.fastq.gz --un fastq_unaligned_TE_consensus/${fq[$i]}.fastq.gz -S fastq_aligned_TE_consensus/${fq[$i]}.sam 2>fastq_unaligned_TE_consensus/${fq[$i]}.map.info
# 	}&
# 	if [[ `expr \( $i + 1 \) % 3` -eq 0 ]]; then
# 		wait
# 	fi
# done


# ######### mirdeep2 expression

# mkdir mirdeep2_out
cd mirdeep2_out
# remove_white_space_in_id.pl /data/users/aimei/shangrui/genomes/dmel-all-chromosome-r6.34.fasta > dmel-all-chromosome-r6.34_no_ws.fasta

# mkdir fasta_file
# for i in ${!fq[@]}
# do
# 	#echo $i
# 	fastq2fasta.pl ../fastq_unaligned_TE_consensus/${fq[$i]}.fastq.gz > fasta_file/${fq[$i]}.fa
# done

nohup mapper.pl config.txt -d -c -i -j -l 18 -m -o 8 -p ~/genomes/Drosophila_melanogaster/small_RNA/genomes/dmel_genome -s reads_collapsed.fa -t reads_collapsed_vs_genome.arf &

nohup quantifier.pl -p ~/genomes/Drosophila_melanogaster/small_RNA/ref_index/dme_hairpin.fa -m ~/genomes/Drosophila_melanogaster/small_RNA/ref_index/dme_miRs.fa -r reads_collapsed.fa -y now -t dme -W &




