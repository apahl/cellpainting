#!/bin/bash
# put just the relevant result files into a tar file
# usage: tar_results <folder (without the "_output" part)>

BASEFOLDER=${1%/}  # removes trailing slash, if there is one
FOLDER=${BASEFOLDER}_output
TARNAME=${BASEFOLDER}_results.tgz

if [[ $BASEFOLDER == "" ]]; then
  echo "Missing parameter."
  echo "usage: tar_results <folder (without the "_output" part)>"
  exit 1
fi

tar cvzf $TARNAME $FOLDER/*.txt $FOLDER/*.tsv $FOLDER/*.cppipe $FOLDER/images
