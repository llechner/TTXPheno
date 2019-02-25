#!/bin/sh

start=1337527
stop=1337553

for i in $(seq $start $stop)
do
    scancel $i
done
