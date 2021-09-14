#!/bin/bash

# check argument number
pcount=$#
if((pcount==0)); then
  echo no args;
  exit;
fi

dir_name=$1
echo directory="${dir_name}"

start=$2
echo start random seed="${start}"

end=$3
echo end random seed="${end}"

config_path=$4
echo config_path="${config_path}"

for((i = ${start}; i < ${end}; i++))
do
  ./MCGRPTW --dir="$dir_name" --random_seed=$i --config="$config_path"
done
