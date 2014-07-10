#!/bin/bash
#set -eu # makes your program exit on error or unbound variable
fldr="/homecentral/srao/Documents/cnrs/simResults/"
spkfile="spkTimes.csv"
fn="spkTimes_theta" 
extn=".csv"
#make

for n in {18..72..36}
do
    echo $n
    ./mysolver $n
done


