#!/bin/bash

SRC="assignment2.cpp"
EXE="assignment2"
PASSES=5
TIMEOUT=10800

DATASETS=("wiki-Vote.txt" "email-Enron.txt" "as-skitter.txt")

g++ -O3 -march=native -std=c++17 -Wall -Wextra "$SRC" -o "$EXE" || exit 1

for DATASET in "${DATASETS[@]}"
do
    OUT_FILE="output_$(basename $DATASET .txt).txt"
    timeout $TIMEOUT ./$EXE "$DATASET" "$PASSES" 2>&1 | tee "$OUT_FILE"
    if [ $? -eq 124 ]; then
        echo "$DATASET NC" | tee -a "$OUT_FILE"
    fi
done
