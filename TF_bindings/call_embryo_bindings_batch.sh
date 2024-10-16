# bed=(`ls Embryo_IDR`)
# for i in ${!bed[@]}
# do
# 	{
# 		perl TF_bindings.pl Embryo_IDR/${bed[$i]} SupplementaryTable1-EmbryoChIP.tsv /public1/users/aimei/genomes/Drosophila_melanogaster/BDGP6_r35/protein_coding.gtf Embryo_bindings
# 	}&
# 	if [[ `expr \( $i + 1 \) % 24` -eq 0 ]]; then
# 		wait
# 	fi
# done
# wait

# bed=(`ls Embryo_IDR | grep -P "CG13296|CG3838|CG7372"`)
# for i in ${!bed[@]}
# do
# 	{
# 		perl TF_bindings.pl Embryo_IDR/${bed[$i]} SupplementaryTable1-EmbryoChIP.tsv /public1/users/aimei/genomes/Drosophila_melanogaster/BDGP6_r35/protein_coding.gtf Embryo_bindings
# 	}&
# 	if [[ `expr \( $i + 1 \) % 24` -eq 0 ]]; then
# 		wait
# 	fi
# done
# wait

bed=(`ls Embryo_IDR | grep -P "CTPsyn"`)
for i in ${!bed[@]}
do
	{
		perl TF_bindings.pl Embryo_IDR/${bed[$i]} SupplementaryTable1-EmbryoChIP.tsv /public1/users/aimei/genomes/Drosophila_melanogaster/BDGP6_r35/protein_coding.gtf Embryo_bindings
	}&
	if [[ `expr \( $i + 1 \) % 24` -eq 0 ]]; then
		wait
	fi
done
wait
