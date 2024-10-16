bed=(`ls W3L_IDR`)
for i in ${!bed[@]}
do
	{
		perl TF_bindings.pl W3L_IDR/${bed[$i]} SupplementaryTable2-W3LChIP.tsv /public1/users/aimei/genomes/Drosophila_melanogaster/BDGP6_r35/protein_coding.gtf W3L_bindings
	}&
	if [[ `expr \( $i + 1 \) % 24` -eq 0 ]]; then
		wait
	fi
done
wait

