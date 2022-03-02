#!/bin/bash

CHECK=1

if [ $# -eq 3 ]; then
  CONFIG_OLD=$1
  CASE_OLD=$2
  CONFIG_NEW=$1
  CASE_NEW=$3
elif [ $# -eq 4 ]; then
  CONFIG_OLD=$1
  CASE_OLD=$2
  CONFIG_NEW=$3
  CASE_NEW=$4
else
  CHECK=0
fi

echo "Copying from $CONFIG_OLD $CASE_OLD to $CONFIG_NEW $CASE_NEW"

if [ ! -d run/nemo_${CONFIG_OLD}_${CASE_OLD} ]; then
  echo "~!@#%^&* ERROR : run/nemo_${CONFIG_OLD}_${CASE_OLD} does not exist"
  echo " "
  CHECK=0
fi

if [ ! -d output/nemo_${CONFIG_OLD}_${CASE_OLD} ]; then
  echo "~!@#%^&* WARNING : output/nemo_${CONFIG_OLD}_${CASE_OLD} does not exist"
  echo "             >>>>>> postprocessing scripts are not being prepared"
  echo " "
fi

if [ -d run/nemo_${CONFIG_NEW}_${CASE_NEW} ]; then
  echo "~!@#%^&* ERROR : run/nemo_${CONFIG_NEW}_${CASE_NEW} does exist"
  echo "             >>>>>> check new config and case names or remove existing directory"
  echo " "
  CHECK=0
fi

if [ $CHECK -ne 1 ]; then
  echo "  Usage: "
  echo "         copy_sim_env_from_existing.sh <existing_config> <existing_case> <new_config> <new_case>"
  echo "         copy_sim_env_from_existing.sh <config> <existing_case> <new_case>"
  echo " "
  echo "  Example: "
  echo "           copy_sim_env_from_existing.sh eORCA025.L75 GNJ001 GNJ002" 
  echo "           copy_sim_env_from_existing.sh AMUXL12.L75 GNJ001 WED12.L121 TEST01" 
  echo " "
  exit
fi

## 1- Prepare running scripts

mkdir -pv run/nemo_${CONFIG_NEW}_${CASE_NEW}

cp -p run/nemo_${CONFIG_OLD}_${CASE_OLD}/*.f90 run/nemo_${CONFIG_NEW}_${CASE_NEW}/.
cp -p run/nemo_${CONFIG_OLD}_${CASE_OLD}/*.xml run/nemo_${CONFIG_NEW}_${CASE_NEW}/.
cp -p run/nemo_${CONFIG_OLD}_${CASE_OLD}/*.sh run/nemo_${CONFIG_NEW}_${CASE_NEW}/.
rm -f run/nemo_${CONFIG_NEW}_${CASE_NEW}/compress_nemo_[0-9]*sh
cp -p run/nemo_${CONFIG_OLD}_${CASE_OLD}/*GENERIC* run/nemo_${CONFIG_NEW}_${CASE_NEW}/.

if [ ! ${CONFIG_OLD} == ${CONFIG_NEW} ]; then
  mv run/nemo_${CONFIG_NEW}_${CASE_NEW}/namelist_nemo_GENERIC_${CONFIG_OLD} run/nemo_${CONFIG_NEW}_${CASE_NEW}/namelist_nemo_GENERIC_${CONFIG_NEW}
  mv run/nemo_${CONFIG_NEW}_${CASE_NEW}/namelist_ice_nemo_GENERIC_${CONFIG_OLD} run/nemo_${CONFIG_NEW}_${CASE_NEW}/namelist_ice_nemo_GENERIC_${CONFIG_NEW}
fi

for file in run/nemo_${CONFIG_OLD}_${CASE_OLD}/run_nemo*sh
do
  file2="run/nemo_${CONFIG_NEW}_${CASE_NEW}/`basename $file`"
  sed -e "s#${CONFIG_OLD}#${CONFIG_NEW}#g ; s#${CASE_OLD}#${CASE_NEW}#g" $file > $file2
  chmod +x $file2
done

## 2- Prepare postprocessing scripts

mkdir -pv output/nemo_${CONFIG_NEW}_${CASE_NEW}

cp -p output/nemo_${CONFIG_OLD}_${CASE_OLD}/*.f90 output/nemo_${CONFIG_NEW}_${CASE_NEW}/. 2>/dev/null || :
cp -p output/nemo_${CONFIG_OLD}_${CASE_OLD}/*.py output/nemo_${CONFIG_NEW}_${CASE_NEW}/. 2>/dev/null || :

for file in output/nemo_${CONFIG_OLD}_${CASE_OLD}/postpro*.sh
do
  file2="output/nemo_${CONFIG_NEW}_${CASE_NEW}/`basename $file`"
  sed -e "s#${CONFIG_OLD}#${CONFIG_NEW}#g ; s#${CASE_OLD}#${CASE_NEW}#g" $file > $file2
  chmod +x $file2
done

## 3- Prepare restart directory

mkdir -pv restart/nemo_${CONFIG_NEW}_${CASE_NEW}
