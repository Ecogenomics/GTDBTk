#!/bin/bash

DATAPATH='database'
for f in $DATAPATH/*.gz
do
  mkdir --parents database/${f:9:3}/${f:13:3}/${f:16:3}/${f:19:3}/ ; mv $f $_
  filef="$(basename -- $f)"
  echo "$filef database/${f:9:3}/${f:13:3}/${f:16:3}/${f:19:3}/ " >> genome_paths.tsv
done