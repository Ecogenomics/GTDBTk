#!/bin/bash

for f in $1/*.gz
do
  filef="$(basename -- $f)"
  mkdir --parents $1/${filef:0:3}/${filef:4:3}/${filef:7:3}/${filef:10:3}/ ; mv $f $_
  echo "$filef $1/${filef:0:3}/${filef:4:3}/${filef:7:3}/${filef:10:3}/" >> genome_paths.tsv

done