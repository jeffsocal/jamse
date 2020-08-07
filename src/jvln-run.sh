#!/bin/bash

D=190522
M="int"

mkdir $D"_"$M"_json"
mkdir $D"_"$M"_results"

cd $D"_"$M"_json"
# mgf2json -p ../mgf/ -m "$M"topic
jvln
jvlnfdr
jvlnid
jvlngrp
cd ../$D"_"$M"_results"
jvlnreport -p ../$D"_"$M"_json"/ -r "YEAST||YEAS||UPI"