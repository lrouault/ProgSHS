#!/bin/bash

echo "1 temperature 2 Eta"
read type

cd fichier
if [ $type == "1" ]
then
    gnuplot -e "a=0; b=100; load '0_Temp.gnu'; pause -1"
fi

if [ $type == "2" ] 
then
    gnuplot -e "a=0; b=100; load '0_Eta.gnu'; pause -1"
fi
cd ..
         
exit 0
      
