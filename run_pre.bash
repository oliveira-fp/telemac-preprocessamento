#!/bin/bash

tempoi=`date`

FF=gfortran
FLAGS='-I/usr/include/ -L/usr/lib -lnetcdff'
ARQUIVOS=`ls *.F90`

./cria_namelist.bash
#echo "EXECUTANDO: ./utm2deg.bash ; ./zz_utm2deg.bash"
#./utm2deg.bash ; ./zz_utm2deg.bash

echo "COMPILANDO SCRIPTS F90"
for arquivo in $ARQUIVOS
do
	$FF $arquivo -o ${arquivo::-4}.o $FLAGS
done

EXECUTAVEIS=`ls *.o`

for executavel in $EXECUTAVEIS
do
	echo "EXECUTANDO: $executavel"
	./$executavel
done

tempof=`date`

echo $tempoi
echo $tempof

