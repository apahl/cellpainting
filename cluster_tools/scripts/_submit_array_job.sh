#!/bin/bash
set -e

FOLDER=$1

while [ $(qstat -u apahl | grep cellprof | wc -l) -gt 0 ]; do
  # do not start a new cellprofiler job when there is already one running
  sleep 600
done

STATUS=$(qsub cellprof_array_96.sge $FOLDER)
# STATUS="Your Job 9367080 has been started."  # for testing

JOB_ID=$(echo $STATUS | awk '{print $3}' | awk 'BEGIN {FS=":"}{print $1}')

echo $JOB_ID $FOLDER
watch_job $JOB_ID $FOLDER

