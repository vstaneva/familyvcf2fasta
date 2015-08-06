#!/bin/bash
set -o nounset
set -o errexit

function pause(){
read -n1 -r -p "$*"
}

FLAGS='--db-attach=yes --leak-check=full --show-reachable=yes'

TOOL="valgrind $FLAGS"
#TOOL="gdb --args"
#TOOL=""

UTESTS=0
TRIO=1

make clean;
make;
chmod 755 test_*
chmod 755 synthetic_trio

if [ $UTESTS -eq 1 ]; then

  $TOOL ./test_phaser;
  echo ""
  #pause "Press any key"
  #echo ""
  #$TOOL ./test_matrix;

fi


if [ $TRIO -eq 1 ]; then
  $TOOL ./synthetic_trio gen_sample/samp21.fa 10 1
  rm -rf LEN*.data
  echo ""
  pause "Press any key"
  echo ""
fi

