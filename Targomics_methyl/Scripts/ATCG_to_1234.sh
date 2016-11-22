

file=$1


cat $file | paste - - - - | awk '{print $2}' > $file"_numeric.tmp"




sed -i "s/A/1,/g" $file"_numeric.tmp"
sed -i "s/T/2,/g" $file"_numeric.tmp"
sed -i "s/C/3,/g" $file"_numeric.tmp"
sed -i "s/G/4,/g" $file"_numeric.tmp"



