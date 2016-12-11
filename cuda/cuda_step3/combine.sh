#!/bin/bash

if [ $# -lt "2" ]; then
    echo "Usage: ./combine.sh <output> <file> [<files>]"
    exit 1
fi

OUTPUT=$1

if [ -f $OUTPUT ]; then
    echo "$OUTPUT is backed up"
    cp $OUTPUT $OUTPUT.bak
fi

shift

echo "" > $OUTPUT

for F in $@; do
    cat $F >> $OUTPUT
    echo "" >> $OUTPUT
done