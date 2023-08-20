#!/bin/bash

echo "Correção do OUTPut" 
for mol in `ls *.log`
	do
       contador=$((contador+1))
       echo "Script no.: " $contador " - Estrutura: " $mol
       input=`echo $mol | sed 's/.log//g'`
       sed '/Dipole orientation:/d' < $input".log" > $input"_temp.log"
       mv $input"_temp.log" $input".log"
       	   
   done
