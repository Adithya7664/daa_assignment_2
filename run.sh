#!/bin/bash

echo "Compiling"
g++ -O3 -march=native daa2.cpp -o bc

echo "Running wiki-Vote"
./bc Wiki-Vote.txt > output_wiki.txt

echo "Running email-Enron"
./bc Email-Enron.txt > output_enron.txt

echo "Running as-skitter"
./bc as-skitter.txt > output_skitter.txt

echo "Done!"
