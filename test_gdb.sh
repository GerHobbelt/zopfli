#!/bin/bash

rm -f ../enwik8.mini.zlib ../enwik8.mini.gz

mkdir -p log

gdb --args ./zopfli  -v ../enwik8.mini  2>&1 | tee log/"zop_out.$(date +%s)"

gunzip ../enwik8.mini.gz -c > testme

diff testme ../enwik8.mini

if [ $? -ne 0 ] ; then 
    echo -e "\nREGRESSION ALERT \n"
else
    echo -e "\n\tno regressions ~yet~\n"
fi 
rm -f testme
