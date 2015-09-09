#!/bin/bash
set -o nounset
set -o errexit

python mfcVCFtoFASTA.py  family_utest.config

for FILE_NAME in "child_1.fa" "child_2.fa" "father_1.fa" "father_2.fa" "mother_1.fa" "mother_2.fa"
do
  if diff ./utest/expected_result/child_1.fa  ./utest/result/child_1.fa > /dev/null
  then
    echo "$FILE_NAME is fine"
  else
    echo "$FILE_NAME differ"
    exit 33
  fi
done
