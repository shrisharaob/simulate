#!/bin/bash
set -eu # makes your program exit on error or unbound variable
fldr="/homecentral/srao/Documents/cnrs/simResults/"
spkfile="spkTimes.csv"
fn="spkTimes_theta" 
extn=".csv"

for n in {126..162..18}
do
    echo $n
    ./mysolver $n
#    mv $fldr$spkfile $fldr$fn$n$extn
done

