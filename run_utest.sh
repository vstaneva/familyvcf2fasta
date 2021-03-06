#!/bin/bash
set -o nounset
set -o errexit

for ID in "1" "2" "3" "4" "5"
do
  echo "***********"
  echo "Runing test n$ID:"
  echo "***********"
  FOLDER="./utest/case_${ID}"
  rm -f ${FOLDER}/result/*
  # N_PATHS = 2 by default
  python mfcVCFtoFASTA.py  ./utest/family_utest_${ID}.config 
  #python mfcVCFtoFASTA.py  family_utest_${ID}.config  1
  mv "phase_string.txt" ${FOLDER}/result/
  for FILE_NAME in "father_1.fa" "father_2.fa" "mother_1.fa" "mother_2.fa" "child_1.fa" "child_2.fa" "phase_string.txt"
  do
    FULL_NAME_RES=${FOLDER}/result/${FILE_NAME}
    FULL_NAME_EXP=${FOLDER}/expected_result/${FILE_NAME}
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
  for FILE_NAME in "child.vcf.new.vcf"
  do
    FULL_NAME_RES=${FOLDER}/${FILE_NAME}
    FULL_NAME_EXP=${FOLDER}/expected_result/${FILE_NAME}
    if diff ${FULL_NAME_RES} ${FULL_NAME_EXP} > /dev/null
    then
      echo "$FILE_NAME equals expected"
    else
      echo "$FILE_NAME differ from expected."
      echo "Please compare"
      echo "${FULL_NAME_RES}" 
      echo "and"
      echo "${FULL_NAME_EXP}" 

      #vimdiff ${FULL_NAME_RES} ${FULL_NAME_EXP}
      exit 33
    fi
  done
  
  echo "test n$ID PASS"
  echo "***********"
done

#LIBERAL MODE:
for ID in "5"
do
  echo "***********"
  echo "Runing test n$ID:"
  echo "***********"
  FOLDER="./utest/case_${ID}"
  rm -f ${FOLDER}/result/*
  # N_PATHS = 2 by default
  #python mfcVCFtoFASTA.py  family_utest_${ID}.config 
  python mfcVCFtoFASTA.py  ./utest/family_utest_${ID}.config  1
  for FILE_NAME in "child.vcf.new.vcf"
  do
    FULL_NAME_RES=${FOLDER}/${FILE_NAME}
    FULL_NAME_EXP=${FOLDER}/expected_result_n1/${FILE_NAME}
    if diff ${FULL_NAME_RES} ${FULL_NAME_EXP} > /dev/null
    then
      echo "$FILE_NAME equals expected"
    else
      echo "$FILE_NAME differ from expected."
      echo "Please compare"
      echo "${FULL_NAME_RES}" 
      echo "and"
      echo "${FULL_NAME_EXP}" 

      #vimdiff ${FULL_NAME_RES} ${FULL_NAME_EXP}
      exit 33
    fi
  done
  
  echo "test n$ID PASS"
  echo "***********"
done

rm -f phase_string.txt 

echo " "
echo "***************************"
echo "Unit tests succesfully run!"
echo "***************************"
