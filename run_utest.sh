#!/bin/bash
set -o nounset
set -o errexit


for ID in "1" 
do
  echo "Test n$ID:"
  echo "***********"
  #python mfcVCFtoFASTA.py  family_utest_1.config
  for FILE_NAME in "child_1.fa" "child_2.fa" "father_1.fa" "father_2.fa" "mother_1.fa" "mother_2.fa"
  do
    if diff ./utest/case_${ID}/expected_result/child_1.fa  ./utest/case_${ID}/result/child_1.fa > /dev/null
    then
      echo "$FILE_NAME equals expected"
    else
      echo "$FILE_NAME differ from expected"
      #vimdiff ./utest/case_${ID}/expected_result/child_1.fa  ./utest/case_${ID}/result/child_1.fa
      exit 33
    fi
  done
done
echo "***************************"
echo "Unit tests succesfully run!"
echo "***************************"
