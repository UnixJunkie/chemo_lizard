#!/bin/bash

rm -f data/test.csv
./molenc.py data/test.smi data/test.csv
cat data/test.csv
