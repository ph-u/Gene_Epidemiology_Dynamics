#!/bin/env bash
# author: ph-u
# script: resMigrate.sh
# desc: copy result files into res/ for download
# in: bash resMigrate.sh
# out: NA
# arg: 0
# date: 20260909

find ../data/ | grep -e "thres\|all\|last" > res.txt

while read L;do
  i=`echo -e ${L} | cut -f 3,6 -d "/" | sed -e "s/\//-/g"`
  cp ${L} ../res/${i}
done < res.txt
rm res.txt

exit
