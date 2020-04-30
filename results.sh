#!/bin/bash
# Call one Zopfli test run
#Num_iterations = $1

# remove previous test
rm -f ../enwik8.mini.zlib ../enwik8.mini.gz

mkdir -p results

# gather all the results in a collection of similarly namedfiles
out_file=results/serial_zopfli_supercloud_results

# for my sanity:
date | tee -a $out_file
# V for verbose
#./zopfli  -v ../enwik8.mini  2>&1 | tee log/"zop_out.$(date +%s)"
./zopfli  --i$1 ../enwik8.mini  2>&1 | tee -a $out_file

# record the size of the output file
ls -l ../enwik8.mini.gz | tee -a $out_file

gunzip ../enwik8.mini.gz -c > testme

diff testme ../enwik8.mini

if [ $? -ne 0 ] ; then 
    echo -e "REGRESSION ALERT\n" | tee -a $out_file
else
    echo -e "Decompression Successful\n" | tee -a $out_file
    rm -f testme
fi 
