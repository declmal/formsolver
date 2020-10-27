#!/bin/bash

function run_unit_tests(){ 
 for file in `ls $1` 
 do 
  ./$1"/"$file 
  OP_MODE=$?
  if [ $OP_MODE -ne 0 ]
  then
    echo $OP_MODE
    exit 1
  fi
  if [ -d $1"/"$file ] 
  then 
   run_unit_tests $1"/"$file 
  fi 
 done 
} 
  
run_unit_tests bin
