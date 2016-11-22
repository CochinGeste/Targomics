
file=$1

sed -i -r "s/NA/F/g" $file


sed -i -r 's/34/R/g' $file
sed -i -r 's/24/L/g' $file
sed -i -r 's/31/X/g' $file

sed -i -r 's/1+/1,/g' $file
sed -i -r 's/2+/2,/g' $file
sed -i -r 's/3+/3,/g' $file
sed -i -r 's/4+/4,/g' $file

sed -i -r 's/R/3,4,/g' $file
sed -i -r 's/L/2,4,/g' $file
sed -i -r 's/X/3,1,/g' $file


sed -i -r "s/F/NA,/g" $file






