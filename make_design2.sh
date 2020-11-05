echo -e 'Sample\tgroup'
for i in "C" "F1" "H" "S0" "S14";
do
	for j in `seq 3`;
		do
		sample="${i}_${j}"
		echo -e ${sample}"\t"${i}
		done
done
