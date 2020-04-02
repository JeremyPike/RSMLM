#!/usr/bin/env bash

OUT_FILE="mem_occ.out"

rm -f $OUT_FILE
touch $OUT_FILE
while :
do
	echo $(date +%H:%M:%S) $(pmap $(pgrep -f test.R) | tail -n 1 | awk '/[0-9]K/{print $2}') >> $OUT_FILE
	sleep 1
done
