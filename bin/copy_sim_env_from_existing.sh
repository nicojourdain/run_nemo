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

if [ ! -d run/nemo_${CONFIG_NEW}_${CASE_NEW} ]; then
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

mkdir run/nemo_${CONFIG_NEW}_${CASE_NEW}

for file in 
calculate_end_date.f90 \\
calculate_end_date_month.f90 \\
compress_nemo_GENERIC.sh \\
domain_def.xml \\
field_def.xml \\
iodef_daily_monthly.xml \\
iodef_daily.xml \\
iodef.xml \\
iodef_monthly_daily.xml \\
namelist_ice_nemo_GENERIC \\
namelist_nemo_GENERIC \\
rebuild.sh \\
run_nemo.sh
do

  if [ -f ${TEMP_NEMO_DIR}/template_run/${file}_${MY_CONFIG} ]; then
    if [ $file == 'run_nemo.sh' ]; then
      sed -e "s#<config>#${MY_CONFIG}#g ; s#<case>#${MY_CASE}#g" ${TEMP_NEMO_DIR}/template_run/${file}_${MY_CONFIG} > run/nemo_${CONFIG}_${CASE}/${file}
    else
      cp -p ${TEMP_NEMO_DIR}/template_run/${file}_${MY_CONFIG} run/nemo_${CONFIG}_${CASE}
    fi
  else
    if [ $file == 'run_nemo.sh' ]; then
      sed -e "s#<config>#${MY_CONFIG}#g ; s#<case>#${MY_CASE}#g" ${TEMP_NEMO_DIR}/template_run/${file} > run/nemo_${CONFIG}_${CASE}/${file}
    else
      cp -p ${TEMP_NEMO_DIR}/template_run/${file} run/nemo_${CONFIG}_${CASE}/.
    fi
  fi

  mv run/nemo_${CONFIG}_${CASE}/namelist_nemo_GENERIC run/nemo_${CONFIG}_${CASE}/namelist_nemo_GENERIC_${MY_CONFIG}
  mv run/nemo_${CONFIG}_${CASE}/namelist_ice_nemo_GENERIC run/nemo_${CONFIG}_${CASE}/namelist_ice_nemo_GENERIC_${MY_CONFIG}

done

## 2- Prepare postprocessing scripts

mkdir output/nemo_${CONFIG}_${CASE}

if [ -f ${TEMP_NEMO_DIR}/template_post/postprocess_nemo_with_monthly_io.sh_${CONFIG} ]; then
  sed -e "s#<config>#${MY_CONFIG}#g ; s#<case>#${MY_CASE}#g" > ${TEMP_NEMO_DIR}/template_post/postprocess_nemo_with_monthly_io.sh_${CONFIG} output/nemo_${CONFIG}_${CASE}/postprocess_nemo_with_monthly_io.sh
else
  sed -e "s#<config>#${MY_CONFIG}#g ; s#<case>#${MY_CASE}#g" > ${TEMP_NEMO_DIR}/template_post/postprocess_nemo_with_monthly_io.sh output/nemo_${CONFIG}_${CASE}/postprocess_nemo_with_monthly_io.sh
fi
 
if [ -f ${TEMP_NEMO_DIR}/template_post/postprocess_nemo.sh_${CONFIG} ]; then
  sed -e "s#<config>#${MY_CONFIG}#g ; s#<case>#${MY_CASE}#g" > ${TEMP_NEMO_DIR}/template_post/postprocess_nemo.sh_${CONFIG} output/nemo_${CONFIG}_${CASE}/postprocess_nemo.sh
else
  sed -e "s#<config>#${MY_CONFIG}#g ; s#<case>#${MY_CASE}#g" > ${TEMP_NEMO_DIR}/template_post/postprocess_nemo.sh output/nemo_${CONFIG}_${CASE}/postprocess_nemo.sh
fi
 
