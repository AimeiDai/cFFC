# 1. TF bindings
cd /public1/users/aimei/cFF_revised/1.TF_bindings
# downloading data
perl downloadList.pl W3LChIPList.txt W3LChIP_2024_4_29_2h_17m.tsv W3L_IDR W3L.narrowbed.downloadlist.sh
perl downloadList.pl EmbryoChIPList.txt Embryo_2024_4_29_4h_38m.tsv Embryo_IDR Embryo.narrowbed.downloadlist.sh

# generating supplementary tables for ChIP-seqs, for those have no FBgnIDs are filled manually
perl generate_tables.pl W3LChIP_2024_4_29_2h_17m.tsv ~/genomes/Drosophila_melanogaster/BDGP6_r35/Drosophila_melanogaster.BDGP6.32.53.gtf SupplementaryTable2-W3LChIP.tsv
perl generate_tables.pl Embryo_2024_4_29_4h_38m.tsv ~/genomes/Drosophila_melanogaster/BDGP6_r35/Drosophila_melanogaster.BDGP6.32.53.gtf SupplementaryTable1-EmbryoChIP.tsv

# identifying targets, those with peaks located within -1.5kb and +500bp of the gene body were identified as TF targets.
# perl TF_bindings.pl Embryo_IDR SupplementaryTable1-EmbryoChIP.tsv ~/genomes/Drosophila_melanogaster/BDGP6_r35/Drosophila_melanogaster.BDGP6.32.53.gtf Embryo_bindings
# perl TF_bindings.pl W3L_IDR SupplementaryTable2-W3LChIP.tsv ~/genomes/Drosophila_melanogaster/BDGP6_r35/Drosophila_melanogaster.BDGP6.32.53.gtf W3L_bindings
cd /public1/users/aimei/cFF_revised/1.TF_bindings
nohup bash call_embryo_bindings_batch.sh &
nohup bash call_W3L_bindings_batch.sh &

cat Embryo_bindings/* | sort | uniq > Embryo_tf_bindings.txt
cat W3L_bindings/* | sort | uniq > W3L_tf_bindings.txt

# 2. miRNA and TF expression at embryo and W3L
cd /public1/users/aimei/cFF_revised/2.miR_TF_expression


# 3. TF RNAi expression
cd /public1/users/aimei/cFF_revised/2.miR_TF_expression/TF_expression
perl  merge_technique_replicates.pl fastq_file /public2/users/aimei/TF_RNAi merging_techrep.sh
perl generate_sample_table_4_DEG.pl sample.list metadata.tsv experiment_report_2024_4_11_9h_55m.tsv sample_table.txt
./hisat2.sh


###### 4. motifs
cd /public1/users/aimei/cFF_revised/4.motifs
# embryonic TFs, conserved miRs, expressed with low threshold
perl find_motifs.pl ../1.TF_bindings/Embryo_tf_bindings.txt Conserved_Family_Info.txt ../2.miR_TF_expression/TF_expression/filtered_embryonic_TFs.txt ../2.miR_TF_expression/miR_expression/embryo_miRs.cons.low.txt ../2.miR_TF_expression/miR_expression/miR_Family_Info_modified.txt embryonic_motifs.consmiR.lowThreshold

# embryonic TFs, conserved miRs, expressed with high threshold
perl find_motifs.pl ../1.TF_bindings/Embryo_tf_bindings.txt Conserved_Family_Info.txt ../2.miR_TF_expression/TF_expression/filtered_embryonic_TFs.txt ../2.miR_TF_expression/miR_expression/embryo_miRs.cons.high.txt ../2.miR_TF_expression/miR_expression/miR_Family_Info_modified.txt embryonic_motifs.consmiR.highThreshold

# embryonic TFs, nonconserved miRs, expressed with low threshold
perl find_motifs.pl ../1.TF_bindings/Embryo_tf_bindings.txt Nonconserved_Family_Info.txt ../2.miR_TF_expression/TF_expression/filtered_embryonic_TFs.txt ../2.miR_TF_expression/miR_expression/embryo_miRs.noncons.low.txt ../2.miR_TF_expression/miR_expression/miR_Family_Info_modified.txt embryonic_motifs.nonconsmiR.lowThreshold

# embryonic TFs, nonconserved miRs, expressed with high threshold
perl find_motifs.pl ../1.TF_bindings/Embryo_tf_bindings.txt Nonconserved_Family_Info.txt ../2.miR_TF_expression/TF_expression/filtered_embryonic_TFs.txt ../2.miR_TF_expression/miR_expression/embryo_miRs.noncons.high.txt ../2.miR_TF_expression/miR_expression/miR_Family_Info_modified.txt embryonic_motifs.nonconsmiR.highThreshold



# L3 TFs, conserved miRs, expressed with low threshold

perl find_motifs.pl ../1.TF_bindings/W3L_tf_bindings.txt Conserved_Family_Info.txt ../2.miR_TF_expression/TF_expression/filtered_L3_TFs.txt ../2.miR_TF_expression/miR_expression/l3_miRs.cons.low.txt ../2.miR_TF_expression/miR_expression/miR_Family_Info_modified.txt l3_motifs.consmiR.lowThreshold

# L3 TFs, conserved miRs, expressed with high threshold
perl find_motifs.pl ../1.TF_bindings/W3L_tf_bindings.txt Conserved_Family_Info.txt ../2.miR_TF_expression/TF_expression/filtered_L3_TFs.txt ../2.miR_TF_expression/miR_expression/l3_miRs.cons.high.txt ../2.miR_TF_expression/miR_expression/miR_Family_Info_modified.txt l3_motifs.consmiR.highThreshold

# L3 TFs, nonconserved miRs, expressed with low threshold
perl find_motifs.pl ../1.TF_bindings/W3L_tf_bindings.txt Nonconserved_Family_Info.txt ../2.miR_TF_expression/TF_expression/filtered_L3_TFs.txt ../2.miR_TF_expression/miR_expression/l3_miRs.noncons.low.txt ../2.miR_TF_expression/miR_expression/miR_Family_Info_modified.txt l3_motifs.nonconsmiR.lowThreshold

# L3 TFs, nonconserved miRs, expressed with high threshold
perl find_motifs.pl ../1.TF_bindings/W3L_tf_bindings.txt Nonconserved_Family_Info.txt ../2.miR_TF_expression/TF_expression/filtered_L3_TFs.txt ../2.miR_TF_expression/miR_expression/l3_miRs.noncons.high.txt ../2.miR_TF_expression/miR_expression/miR_Family_Info_modified.txt l3_motifs.nonconsmiR.highThreshold


######### redundancies and conservation
cd /public1/users/aimei/cFF_revised/4.motifs
embryo_consmiR.lowThreshold.redundancy.R

####### polymorphisms and selection
cd /public1/users/aimei/cFF_revised/5.polymorphism
# convert target sites coordinates into genomic coordinate
# W3L, conserved miR, low threshold
perl convert_UTR_coordinate.pl ../4.motifs/l3_consmiR.lowThreshold.redundancy TSFly_r6.19_3UTRs.gff l3_consmiR.lowThreshold.UTRcoor
msort -k d1 -k n2 l3_consmiR.lowThreshold.UTRcoor | uniq > l3_consmiR.lowThreshold.UTRcoor_sorted
#divergence
perl 2.1.count_site_divergence.pl wholeGenomeDivPolDmel/dme.whole.div l3_consmiR.lowThreshold.UTRcoor_sorted >> l3_consmiR.lowThreshold.mkt_ds.out_by_AG &
#polymorphism
perl 2.2.count_site_polymorphism.pl l3_consmiR.lowThreshold.UTRcoor_sorted wholeGenomeDivPolDmel/dme.whole.pol 30 l3_consmiR.lowThreshold.redundancy_DAF.tab >> l3_consmiR.lowThreshold.mkt_ds.out_by_AG &

# embryo, conserved miR, low threshold
perl convert_UTR_coordinate.pl ../4.motifs/Embryonic_consmiR.lowThreshold.redundancy TSFly_r6.19_3UTRs.gff Embryonic_consmiR.lowThreshold.UTRcoor
msort -k d1 -k n2 Embryonic_consmiR.lowThreshold.UTRcoor | uniq > Embryonic_consmiR.lowThreshold.UTRcoor_sorted
#divergence
perl 2.1.count_site_divergence.pl wholeGenomeDivPolDmel/dme.whole.div Embryonic_consmiR.lowThreshold.UTRcoor_sorted >> Embryonic_consmiR.lowThreshold.mkt_ds.out_by_AG &
#polymorphism
perl 2.2.count_site_polymorphism.pl Embryonic_consmiR.lowThreshold.UTRcoor_sorted wholeGenomeDivPolDmel/dme.whole.pol 30 Embryonic_consmiR.lowThreshold.redundancy_DAF.tab >> Embryonic_consmiR.lowThreshold.mkt_ds.out_by_AG &

#############################################
# W3L, conserved miR, high threshold
perl convert_UTR_coordinate.pl ../4.motifs/l3_consmiR.highThreshold.redundancy TSFly_r6.19_3UTRs.gff l3_consmiR.highThreshold.UTRcoor
msort -k d1 -k n2 l3_consmiR.highThreshold.UTRcoor | uniq > l3_consmiR.highThreshold.UTRcoor_sorted
#divergence
perl 2.1.count_site_divergence.pl wholeGenomeDivPolDmel/dme.whole.div l3_consmiR.highThreshold.UTRcoor_sorted >> l3_consmiR.highThreshold.mkt_ds.out_by_AG &
#polymorphism
perl 2.2.count_site_polymorphism.pl l3_consmiR.highThreshold.UTRcoor_sorted wholeGenomeDivPolDmel/dme.whole.pol 30 l3_consmiR.highThreshold.redundancy_DAF.tab >> l3_consmiR.highThreshold.mkt_ds.out_by_AG &

# embryo, conserved miR, high threshold
perl convert_UTR_coordinate.pl ../4.motifs/Embryonic_consmiR.highThreshold.redundancy TSFly_r6.19_3UTRs.gff Embryonic_consmiR.highThreshold.UTRcoor
msort -k d1 -k n2 Embryonic_consmiR.highThreshold.UTRcoor | uniq > Embryonic_consmiR.highThreshold.UTRcoor_sorted
#divergence
perl 2.1.count_site_divergence.pl wholeGenomeDivPolDmel/dme.whole.div Embryonic_consmiR.highThreshold.UTRcoor_sorted >> Embryonic_consmiR.highThreshold.mkt_ds.out_by_AG &
#polymorphism
perl 2.2.count_site_polymorphism.pl Embryonic_consmiR.highThreshold.UTRcoor_sorted wholeGenomeDivPolDmel/dme.whole.pol 30 Embryonic_consmiR.highThreshold.redundancy_DAF.tab >> Embryonic_consmiR.highThreshold.mkt_ds.out_by_AG &


# W3L, nonconserved miR, low threshold
perl convert_UTR_coordinate.pl ../4.motifs/l3_nonconsmiR.lowThreshold.redundancy TSFly_r6.19_3UTRs.gff l3_nonconsmiR.lowThreshold.UTRcoor
msort -k d1 -k n2 l3_nonconsmiR.lowThreshold.UTRcoor | uniq > l3_nonconsmiR.lowThreshold.UTRcoor_sorted
#divergence
perl 2.1.count_site_divergence.pl wholeGenomeDivPolDmel/dme.whole.div l3_nonconsmiR.lowThreshold.UTRcoor_sorted >> l3_nonconsmiR.lowThreshold.mkt_ds.out_by_AG &
#polymorphism
perl 2.2.count_site_polymorphism.pl l3_nonconsmiR.lowThreshold.UTRcoor_sorted wholeGenomeDivPolDmel/dme.whole.pol 30 l3_nonconsmiR.lowThreshold.redundancy_DAF.tab >> l3_nonconsmiR.lowThreshold.mkt_ds.out_by_AG &

# embryo, nonconserved miR, low threshold
perl convert_UTR_coordinate.pl ../4.motifs/Embryonic_nonconsmiR.lowThreshold.redundancy TSFly_r6.19_3UTRs.gff Embryonic_nonconsmiR.lowThreshold.UTRcoor
msort -k d1 -k n2 Embryonic_nonconsmiR.lowThreshold.UTRcoor | uniq > Embryonic_nonconsmiR.lowThreshold.UTRcoor_sorted
#divergence
perl 2.1.count_site_divergence.pl wholeGenomeDivPolDmel/dme.whole.div Embryonic_nonconsmiR.lowThreshold.UTRcoor_sorted >> Embryonic_nonconsmiR.lowThreshold.mkt_ds.out_by_AG &
#polymorphism
perl 2.2.count_site_polymorphism.pl Embryonic_nonconsmiR.lowThreshold.UTRcoor_sorted wholeGenomeDivPolDmel/dme.whole.pol 30 Embryonic_nonconsmiR.lowThreshold.redundancy_DAF.tab >> Embryonic_nonconsmiR.lowThreshold.mkt_ds.out_by_AG &
