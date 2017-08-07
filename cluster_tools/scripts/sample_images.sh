#!/bin/bash
# usage: sample_images PLATE
set -e

PLATE=${1%/}
SRC_DIR=/ptmp/allg/apahl/cp/queue/$PLATE
DEST_DIR=/ptmp/allg/apahl/cp/output/${PLATE}_output

if [ -z "$SRC_DIR" ] || [ ! -d "$SRC_DIR" ]; then
  echo Source Dir $SRC_DIR does not exist
  exit 1
fi

if [ -z $DEST_DIR ] || [ ! -d $DEST_DIR ]; then
  echo Destination Dir $DEST_DIR does not exist
  exit 1
fi

mkdir -p $DEST_DIR/images
for row in {A..P}; do
  for col_no in {1..24}; do
    col=$(printf "%02d" $col_no)
    if [ $col -ne 11 -a $col -ne 12 ] || [ $row$col == H11 ]; then
      for w in {1..5}; do
        for f in "$SRC_DIR"/*/*/*$row${col}_s5_w$w*.tif; do
          if [ -e "$f" ] && [[ "$f" != *_thumb* ]]; then
            echo converting $PLATE : $row${col} channel $w ...
            set +e
            convert "$f" -resize 750x750 \
                   -brightness-contrast 60x70 \
                   $DEST_DIR/images/$row${col}_w$w.jpg > /dev/null 2>&1
            set -e
          fi
        done
      done
    fi
  done
done
