file="../YeastSequences.txt"
window=73

paramset="../cgDNA+_ps1.txt"

module load matlab

count=1;
while IFS= read seq
do
  # echo "$seq"
  seqlen=${#seq}
  # echo "$seqlen" #cia komentaras
  
  n=$(expr $seqlen - $window + 1)

  for ((i = 1; i <= $n; i=i+25)); do
    rm my_seq.txt
    rm my_seq_data.txt
    echo ${seq:i-1:window} >> my_seq.txt
    #echo ${seq:i-1:window}
    ../build_seq_data  -p $paramset  -s my_seq.txt  -o my_seq_data.txt
    ../run_cgDNAmc_from_seq  -e t0,s1,r,t1,t11  -l Output  -i my_seq_data.txt  -p $paramset  -a 100000  -d 10  -j n
    matlab -nosplash -nodisplay -nodesktop -r "run('/home/ginisnai/cgDNApmc/cgDNAplus_mc_parallel/cgDNApmc/runmatlab($count)'); exit; "
done
count=$count+1;
done <"$file"

