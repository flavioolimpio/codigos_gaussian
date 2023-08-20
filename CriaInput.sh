for z in PROD1A PROD2 PROD3 PROD4A PROD4B PROD5 PROD6 PROD7 PROD8 PROD9 TS_PROD1 TS_PROD2 TS_PROD3 TS_PROD4 TS_PROD5 TS_PROD6 TS_PROD7 TS_PROD8 TS_PROD9; do
for x in M06HF; do
for y in 6-31g 6-31ga 6-31mg 6-31mga  6-311g 6-311mg 6-311mga 6-311mmgaa cc-pvdz; do
    if [ "${y}" = "6-311mmgaa" ]; then
          y1="6-311++g**"
        elif [ "${y}" = "6-31mga" ]; then
           y1="6-31+G*"
		elif [ "${y}" = "6-31mg" ]; then
            y1="6-31+G"	
		elif [ "${y}" = "6-311mg" ]; then
            y1="6-311+G"	
        elif [ "${y}" = "6-31ga" ]; then
            y1="6-31G*"
        elif [ "${y}" = "6-311mga" ]; then
           y1="6-311+G*"
    else y1=${y}
    fi 
		input=${z}"_"${x}"_"${y} 
		cp ${z}".chk" ${input}".chk"
        echo "%nprocshared=4" >> ${input}".gjf"
        echo "%mem=8GB"      >> ${input}".gjf"
        echo "%chk="${input}".chk" >> ${input}".gjf" 
		if [ "${z}" = "PROD1A" ] || [ "${z}" = "PROD2" ] || [ "${z}" = "PROD3" ] || [ "${z}" = "PROD4A" ] || [ "${z}" = "PROD4B" ] || [ "${z}" = "PROD5" ] || [ "${z}" = "PROD6" ] || [ "${z}" = "PROD7" ] || [ "${z}" = "PROD8" ] || [ "${z}" = "PROD9" ]; then
			echo "# opt=(readfc,maxcycle=100000,noeigentest) guess=read int=ultrafine freq=noraman ${x}/${y1} geom=(allcheck) scf=(novaracc,xqc,maxcycle=10000) scrf=smd empiricaldispersion(pfd)"  >> ${input}".gjf"
		elif [ "${z}" = "PROD1A" ] || [ "${z}" = "PROD2" ] || [ "${z}" = "PROD3" ] || [ "${z}" = "PROD4A" ] || [ "${z}" = "PROD4B" ] || [ "${z}" = "PROD5" ] || [ "${z}" = "TS_PROD6" ] || [ "${z}" = "TS_PROD7" ] || [ "${z}" = "TS_PROD8" ] || [ "${z}" = "TS_PROD9" ]; then
			echo "# opt=(readfc,modredundant,maxcycle=100000,noeigentest) guess=read int=ultrafine freq=noraman ${x}/${y1} geom=(allcheck) scf=(novaracc,xqc,maxcycle=10000) scrf=smd empiricaldispersion(pfd)"  >> ${input}".gjf"
		fi
		echo " " >> ${input}".gjf"
        cat ${z}".dat" >> ${input}".gjf"
        echo " " >> ${input}".gjf"
		echo " " >> ${input}".gjf"
		echo "g16 ${z}"_"${x}"_"${y}.gjf" >> "fila_${z}.sh"			
done
done
done

chmod 755 fila* 

echo "inputs finalizados!!!"
