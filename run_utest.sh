#!/bin/bash
set -o nounset
set -o errexit


for ID in "1"  "2"
do
  echo "***********"
  echo "Runing test n$ID:"
  echo "***********"
  python mfcVCFtoFASTA.py  family_utest_${ID}.config
  FOLDER="./utest/case_${ID}"
  for FILE_NAME in "child_1.fa" "child_2.fa" "father_1.fa" "father_2.fa" "mother_1.fa" "mother_2.fa"
  do
    FULL_NAME_RES=${FOLDER}/result/${FILE_NAME}
    FULL_NAME_EXP=${FOLDER}/expected_result/${FILE_NAME}
    #if diff ./utest/case_${ID}/expected_result/child_1.fa  ./utest/case_${ID}/result/child_1.fa > /dev/null
    if diff ${FULL_NAME_RES} ${FULL_NAME_EXP} > /dev/null
    then
      echo "$FILE_NAME equals expected"
    else
      echo "$FILE_NAME differ from expected."
      echo "Please compare"
      echo "${FULL_NAME_RES}" 
      echo "and"
      echo "${FULL_NAME_EXP}" 

      vimdiff ${FULL_NAME_RES} ${FULL_NAME_EXP}
      exit 33
    fi
  done
  echo "test n$ID PASS"
  echo "***********"
done
echo " "
echo "***************************"
echo "Unit tests succesfully run!"
echo "***************************"
