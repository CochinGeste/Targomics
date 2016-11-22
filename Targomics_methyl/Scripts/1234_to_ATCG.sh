

file=$1


sed -i "s/1 /A,/g" $file
sed -i "s/2 /T,/g" $file
sed -i "s/3 /C,/g" $file
sed -i "s/4 /G,/g" $file

sed -i "s/1/A,/g" $file
sed -i "s/2/T,/g" $file
sed -i "s/3/C,/g" $file
sed -i "s/4/G,/g" $file


sed -i "s/NA /NA,/g" $file


