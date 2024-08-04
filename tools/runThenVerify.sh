#!/bin/bash

trap 'echo "Script interrupted"; exit 1' INT

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <input_file> <output_path> [options]"
    exit 1
fi

INPUT_FILE=$1
OUTPUT_PATH=$2
shift 2
OPTIONS=$@

BASENAME=$(basename "$INPUT_FILE" .aig)

mkdir -p "$OUTPUT_PATH"

./miniCAR $OPTIONS "$INPUT_FILE" "$OUTPUT_PATH"

../aiger/aigsim -c "$INPUT_FILE" "$OUTPUT_PATH/$BASENAME.res"

if [ $? -eq 0 ]; then
    FIRST_LINE=$(head -n 1 "$OUTPUT_PATH/$BASENAME.res")
    if [ "$FIRST_LINE" == "1" ]; then
    echo "Verification passed: UNSAFE"
    else
    echo "Verification passed: SAFE"
    fi
else
    echo "Verification failed"
fi
