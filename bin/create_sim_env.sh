#!/bin/bash

if [ $# -ne 2 ]; then
  echo "  Usage: create_sim_env.sh <config> <case>"
  echo " "
  echo "  Example: create_sim_env.sh eORCA025.L75 GNJ001" 
  echo " "
  exit
fi

CONFIG=$1
CASE=$2

## 1- Prepare running scripts

if [ -d run/nemo_${CONFIG}_${CASE} ]; then
  echo "~!@#$%& ERROR: run/nemo_${CONFIG}_${CASE} already exist >>> delete or change config/case names"
  exit
fi
mkdir run/nemo_${CONFIG}_${CASE}

for file in calculate_end_date.f90 \
calculate_end_date_month.f90 \
compress_nemo_GENERIC.sh \
domain_def.xml \
field_def.xml \
iodef_daily_monthly.xml \
iodef_daily.xml \
iodef_monthly_daily.xml \
namelist_ice_nemo_GENERIC \
namelist_nemo_GENERIC \
rebuild.sh \
run_nemo.sh
do

  if [ -f ${TEMP_NEMO_DIR}/template_run/${file}_${CONFIG} ]; then
    if [ $file == 'run_nemo.sh' ]; then
      sed -e "s#<config>#${CONFIG}#g ; s#<case>#${CASE}#g" ${TEMP_NEMO_DIR}/template_run/${file}_${CONFIG} > run/nemo_${CONFIG}_${CASE}/${file}
      chmod +x run/nemo_${CONFIG}_${CASE}/${file}
    else
      cp -p ${TEMP_NEMO_DIR}/template_run/${file}_${CONFIG} run/nemo_${CONFIG}_${CASE}
    fi
    echo "Copying template_run/${file}_${CONFIG}"
  else
    if [ $file == 'run_nemo.sh' ]; then
      sed -e "s#<config>#${CONFIG}#g ; s#<case>#${CASE}#g" ${TEMP_NEMO_DIR}/template_run/${file} > run/nemo_${CONFIG}_${CASE}/${file}
      chmod +x run/nemo_${CONFIG}_${CASE}/${file}
    else
      cp -p ${TEMP_NEMO_DIR}/template_run/${file} run/nemo_${CONFIG}_${CASE}/.
    fi
    echo "Copying template_run/${file} (standard template)"
  fi

done

mv run/nemo_${CONFIG}_${CASE}/namelist_nemo_GENERIC run/nemo_${CONFIG}_${CASE}/namelist_nemo_GENERIC_${CONFIG}
mv run/nemo_${CONFIG}_${CASE}/namelist_ice_nemo_GENERIC run/nemo_${CONFIG}_${CASE}/namelist_ice_nemo_GENERIC_${CONFIG}

## 2- Prepare postprocessing scripts

mkdir output/nemo_${CONFIG}_${CASE}

if [ -f ${TEMP_NEMO_DIR}/template_post/postprocess_nemo_with_monthly_io.sh_${CONFIG} ]; then
  sed -e "s#<config>#${CONFIG}#g ; s#<case>#${CASE}#g" ${TEMP_NEMO_DIR}/template_post/postprocess_nemo_with_monthly_io.sh_${CONFIG} > output/nemo_${CONFIG}_${CASE}/postprocess_nemo_with_monthly_io.sh
  echo "Copying template_post/postprocess_nemo_with_monthly_io.sh_${CONFIG}"
else
  sed -e "s#<config>#${CONFIG}#g ; s#<case>#${CASE}#g" ${TEMP_NEMO_DIR}/template_post/postprocess_nemo_with_monthly_io.sh > output/nemo_${CONFIG}_${CASE}/postprocess_nemo_with_monthly_io.sh
  echo "Copying template_post/postprocess_nemo_with_monthly_io.sh (standard template)"
fi
chmod +x output/nemo_${CONFIG}_${CASE}/postprocess_nemo_with_monthly_io.sh
 
if [ -f ${TEMP_NEMO_DIR}/template_post/postprocess_nemo.sh_${CONFIG} ]; then
  sed -e "s#<config>#${CONFIG}#g ; s#<case>#${CASE}#g" ${TEMP_NEMO_DIR}/template_post/postprocess_nemo.sh_${CONFIG} > output/nemo_${CONFIG}_${CASE}/postprocess_nemo.sh
  echo "Copying template_post/postprocess_nemo.sh_${CONFIG}"
else
  sed -e "s#<config>#${CONFIG}#g ; s#<case>#${CASE}#g" ${TEMP_NEMO_DIR}/template_post/postprocess_nemo.sh > output/nemo_${CONFIG}_${CASE}/postprocess_nemo.sh
  echo "Copying template_post/postprocess_nemo.sh (standard template)"
fi
chmod +x output/nemo_${CONFIG}_${CASE}/postprocess_nemo.sh 

