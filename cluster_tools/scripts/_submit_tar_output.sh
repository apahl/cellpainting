#!/bin/bash
set -e

FOLDER=$1

while [ $(qstat -u apahl | grep tar_out | wc -l) -gt 0 ]; do
  # do not start a new tar job when there is already one running
  sleep 300
done

qsub tar_output.sge $FOLDER
