#!/bin/bash
# run.sh <hailfinder.bif> <records.dat>
if [ $# -lt 2 ]; then
  echo "Usage: $0 <hailfinder.bif> <records.dat>"
  exit 1
fi

BIF="$1"
DATA="$2"

if [ ! -f "$BIF" ]; then
  echo "BIF file not found: $BIF"
  exit 1
fi

if [ ! -f "$DATA" ]; then
  echo "Data file not found: $DATA"
  exit 1
fi

# run the program; output file is solved_hailfinder.bif
./learn_bn "$BIF" "$DATA"
