#!/bin/bash

if [ $# -ne 4 ]; then
    echo "Wrong number of parameters. [progname] template outconfigbasename outlogbasename nbRuns"
    exit
fi
tmplate=$1
outbasename=$2
outlogbasename=$3
nbRuns=$4

nbRob='100'

#sigma='0.1'

selpres='1.0'

broadcast='0 2 -1'
frequence='10'

distance='1000 200 80 40'



listProp=`parallel --header : echo R{1}.SP{2}.FR{3}.BC{4}.D{5} gInitialNumberOfRobots={f1} gSelPressure={f2} gBroadcastTime={f3} gMatingOperator={f4} gCommunicationRange={f5} ::: f1 $nbRob ::: f2 $selpres ::: f3 $frequence ::: f4 $broadcast ::: f5 $distance`

#listProp=`parallel --header : echo R{1}.T{2}.B{3}.S{4} gInitialNumberOfRobots={f1} gTaskSeq={f2} gBrait={f3} gSigmaRef={f4} ::: f1 $nbRob ::: f2 $taskSeq ::: f3 $ctrlSetup ::: f4 $sigma`

#echo "$listProp"
#exit

program="./roborobo -l "

mkdir $outlogbasename
commandFile=$tmplate.parallel

rm $commandFile
touch $commandFile

while read -r line
do
suffix=`echo $line | cut -f1 -d" " | sed -e 's/,-1//' |sed -e 's/,/T/g'`
#echo $suffix 

cp $tmplate $outbasename-$suffix.properties
mkdir $outlogbasename/$suffix
#echo $line
echo $line | perl -p -e  's/^.*? //' | sed -e 's/ /\n/g' >> $outbasename-$suffix.properties



    for (( j=1; j<=$nbRuns; j++))
    do
	echo "$program $outbasename-$suffix.properties 1> $outlogbasename/$suffix/run-$j.log" >> $commandFile    
    done

done <<< "$listProp"


#parallel -j15 -a $commandFile



#./tools/template-launchMultiple.sh config/tColl2/template-medea config/tColl2/medea-colors logs/expColors 30


#for i in * ; do paste $i/*.log > $i.all.log; done

#for i in */ ; do paste $i/*.log.1 > $i/all.log.1; done

#for i in *.log ; do tail -n +21 $i > tmp ; rm $i ; mv tmp $i ; done ;
#for i in * ; do for j in $i/*.log ; do tail -n +22 $j > tmp ; rm $j ; mv tmp $j ; done ; done ;

#for i in *.log ; do head -n -6 $i > tmp ; rm $i ; mv tmp $i ; done ;
#for i in * ; do for j in $i/*.log ; do head -n -6 $j > tmp ; rm $j ; mv tmp $j ; done ; done ; 

#scp -rf -P 3003 cecyle.irit.fr:~/roborobo3/logs/collect2_with_positions_*/ .


# for i in */ ; do for j in $i/*.log ; do ../cutIntoColumnFiles.sh $j ;  done ; done ; 
