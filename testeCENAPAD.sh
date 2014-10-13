#PBS -q par22
#PBS -N histPAR
#PBS -l ncpus=2
#PBS -m abe
#PBS -V
#!/bin/sh
echo "----------------------------------------"
echo "Inicio do job:" `date`
echo "Hostname: " `hostname`
echo "PWD: " $PWD
outputFolder="/home/sgi/proj/proj595/borges/Documents/exomas_teste/density"
inputFolder="/home/sgi/proj/proj595/borges/Documents/exomas_teste/dep"
cd $inputFolder
files=`ls *dep`
for file in $files
do
LC_ALL=C locale
orderDEP=(`cut -f 3 $file | sort -n | uniq`) 
cut -f 3 $file | sort -n  > $file.tmp 
median=`cat $file.tmp | awk '{arr[NR]=$1} END { if (NR%2==1) print arr[(NR+1)/2]; else print (arr[NR/2]+arr[NR/2+1])/2}'` 
sizeOfFileDep=`wc -l $file | awk '{print $1}'` 
sizeOforderDEP=`echo ${#orderDEP[@]}`
totalCounts=`awk '{s+=$1}END{print s}' $file.tmp`
efficiencyCount=$totalCounts
for i in `seq 1 $sizeOforderDEP`
do
let i=$i-1
dep=${orderDEP[$i]}
cont=`grep -w $dep $file.tmp -c`
proportion=`bc -l <<< 100*$cont*$dep/$totalCounts`
if [ $dep = 0 ]
then
proportion=`bc -l <<< 100*$cont/$totalCounts`
fi
efficiencyPorportion=`bc -l <<< 100*$efficiencyCount/$totalCounts`
echo $dep $cont $proportion $efficiencyPorportion
efficiencyCount=`bc -l <<< $efficiencyCount-$cont*$dep`
done > $outputFolder/$file-$median-hist
rm $file.tmp
done
echo "Final do job:" `date`
echo "----------------------------------------"
