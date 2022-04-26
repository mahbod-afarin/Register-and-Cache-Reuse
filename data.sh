export MKL_ENABLE_INSTRUCTIONS=AVX2
export MKL_NUM_THREADS=1

cd build
output_name=$1

./reg_reuse &> $output_name
./cache_part3 >> $output_name

echo "" >> $output_name
echo "*********** Cache reuse - Part 4. ***********" >> $output_name

declare -a opt_flags=("0" "1" "2" "3")
opt_nums=${#opt_flags[@]}

for (( i=1; i<${opt_nums}+1; i++ ));
do
  out_string="O${opt_flags[$i-1]} : \c"
  echo -e $out_string >> $output_name
  executable_name="cache_part4_o${opt_flags[$i-1]}"
  ./$executable_name >> $output_name
done

cd ../ && cp build/$output_name data/$output_name


