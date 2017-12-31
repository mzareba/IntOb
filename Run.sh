#!/bin/bash

rm -rf optimizedEquationNo*

#mkdir ./Results
cp equationList.txt backup_equationList.txt

let linesMax=`cat -n equationList.txt | tail -1 | awk '{print $1}'`

#Main loop
let step=1
for i in `seq 1 $linesMax` ;do
IFS=$'\n'
lines=($(cat equationList.txt))
let line=0

#Write upper part R script
echo ${lines[$line]} >> equation.txt

#Make ready to use R script
name=`echo "optimizedEquationNo_"$step`
mkdir ./$name

#Run R
R CMD BATCH --slave ./farmacja_v1.R

#End process
mv *.Rout ./$name
mv Results.txt ./$name
mv optimizedEquation.txt ./$name
mv equation.txt ./$name
#tar -cvf $name.tar.gz ./$name

#rm -r ./$name
#mv *.tar.gz ./Results
tail -n +2 equationList.txt > temp.txt
cat temp.txt > equationList.txt
rm -r temp.txt
rm -r equation.txt
let step=step+1

done


cp -f backup_equationList.txt equationList.txt
exit 0

