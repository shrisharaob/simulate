#!/bin/bash

clusters={75..83..1} 
minmem=22948084
for kCkuster in clusters, do
    ssh srao@neurphys$kCkuster
    freeMem=`cat /proc/meminfo | grep MemFree | grep "[0-9]{1,100}"`
    if [ freeMem >= minmem ]; then
        screen -dmS sh -c "echo $PATH; exec /bin/sh"
    fi
#    screen -dmS cluster$kCkuster /homecentral/Documents/code/mysolver  
#reen -dmS tst03 sh -c "./mysolver2; exec /bin/bash"
