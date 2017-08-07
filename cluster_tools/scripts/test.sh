#!/bin/bash
for col_no in {1..24}; do
    col=$(printf "%02d" $col_no)
    echo $col
done

