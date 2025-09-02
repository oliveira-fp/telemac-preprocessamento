#!/bin/bash

tempoi=`date`

FF=gfortran
FLAGS='-I/usr/include/ -L/usr/lib -lnetcdff'

ARQUIVOS=`ls 10_*.F90 11_*.F90 12_*.F90 13_*.F90`

echo -e "\nCOMPILANDO SCRIPTS F90\n"
for arquivo in $ARQUIVOS
do
        echo "COMPILANDO: $arquivo"
	$FF $arquivo -o ${arquivo::-4}.o $FLAGS
done

EXECUTAVEIS=`ls *.o`
echo -e "\nEXECUTANDO SCRIPTS F90\n"
for executavel in $EXECUTAVEIS
do
        echo -e "\n#####################################################\n"
        echo -e "\nEXECUTANDO: $executavel\n"
	./$executavel
done
rm *.o

tempof=`date`

echo -e "\n#####################################################\n"
echo 'EXECUÇÃO DO PREPROCESSAMENTO'
echo -e "PARTE 3 - INTERPOLAÇÃO DOS DADOS ATMOSFÉRICOS E OCEÂNICOS\nPARA GRADE DO MODELO TELEMAC"
echo -e "\nINÍCIO: ${tempoi}"
echo -e "FIM   : ${tempof}"
echo -e "\n#####################################################\n"
