#! /bin/bash

# store the current date
date=$(date +%Y-%m-%d)

# download all current proteins from fpbase
python3 \
    gather_fpbase.py \
    > fpbase-proteins_"${date}".tsv

