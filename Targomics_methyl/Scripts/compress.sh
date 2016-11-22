
file=$1


sed -i -r "s/1/A/g" $file
sed -i -r "s/2/T/g" $file
sed -i -r "s/3/C/g" $file
sed -i -r "s/4/G/g" $file



sed -i -r 's/"//g' $file
sed -i -r 's/NA//g' $file

sed -i -r 's/CA/P/g' $file
sed -i -r 's/CG/O/g' $file
sed -i -r 's/TG/L/g' $file

sed -i -r 's/T+/T/g' $file
sed -i -r 's/C+/C/g' $file
sed -i -r 's/A+/A/g' $file
sed -i -r 's/G+/G/g' $file

sed -i -r 's/P/CA/g' $file
sed -i -r 's/O/CG/g' $file
sed -i -r 's/L/TG/g' $file

sed -i -r "s/A/1,/g" $file
sed -i -r "s/T/2,/g" $file
sed -i -r "s/C/3,/g" $file
sed -i -r "s/G/4,/g" $file

