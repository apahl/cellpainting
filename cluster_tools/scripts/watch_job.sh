#!/bin/bash
# send a mail when the watched array job finishes.
# usage: watch_job <job_id> <folder>

JOB_ID=$1
PLATE=${2%/}  # removes trailing slash, if there is one

if [[ $JOB_ID == "" || $PLATE == "" ]]; then
  echo "Missing parameter."
  echo "usage: watch_job <job_id> <folder>"
  exit 1
fi

JOB_LOG=/ptmp/allg/apahl/job_${JOB_ID}_${PLATE}.log

while true; do
  if [ -e $JOB_LOG ]; then
    NUM_FINISHED=$(cat $JOB_LOG | wc -l)
    if [ $NUM_FINISHED -eq 96 ]; then
      NUM_FAILED=$(grep FAILED $JOB_LOG | wc -l)
      if [ $NUM_FAILED -eq 0 ]; then
        STATUS="All tasks finished without error."  
      else
        STATUS="$NUM_FAILED tasks finished with error."
      fi
      echo -e "Job $JOB_ID cellprof_array_96 $PLATE has finished.\n$STATUS" | mail -s "Job finished" axel.pahl@mpi-dortmund.mpg.de
      exit 0
    fi
  fi
  sleep 600
done
