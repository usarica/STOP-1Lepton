#!/bin/sh

INFILE=$1
SCRIPTNAME="doFullProduction.cc"
FCN="doFullProduction"
QUEUE="default"
DATE=$(date +%y%m%d)
OUTDIR="./output/${DATE}"

mkdir -p $OUTDIR

while IFS='' read -r line || [[ -n "$line" ]]; do
  submitSTOP1LGenericJob.sh "$SCRIPTNAME" "$FCN" "$line" "$QUEUE" "$OUTDIR"
done < "$INFILE"
