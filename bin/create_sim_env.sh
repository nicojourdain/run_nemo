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

if [ `hostname -d |grep occigen |wc -c` -gt 0 ]; then
  HOST='occigen'
elif [ `hostname -d |grep irene |wc -c` -gt 0 ]; then
  HOST='irene'
else
  echo '~!@#$%^&* ERROR : you need to define header for this machine'
  echo "`hostname -d` not defined (only occigen and irene available)"
  echo '>>>>>>>>>>>> STOP !!!'
  exit
fi
echo "###################################################"
echo "  Host is $HOST: choosing appropriate headers "
echo "###################################################"

## 1- Prepare running scripts

if [ -d run/nemo_${CONFIG}_${CASE} ]; then
  echo "~!@#$%& ERROR: run/nemo_${CONFIG}_${CASE} already exists >>> delete or change config/case names"
  exit
fi
mkdir run/nemo_${CONFIG}_${CASE}

for file in calculate_end_date.f90 \
calculate_end_date_month.f90 \
compress_nemo_GENERIC.sh \
context_nemo.xml \
domain_def_nemo.xml \
field_def_nemo-ice.xml \
field_def_nemo-oce.xml \
file_def_nemo-ice_1d.xml \
file_def_nemo-ice_1d1m.xml \
file_def_nemo-ice_1m.xml \
file_def_nemo-oce_1d.xml \
file_def_nemo-oce_1d1m.xml \
file_def_nemo-oce_1m.xml \
grid_def_nemo.xml \
iodef.xml \
namelist_nemo-ice_GENERIC \
namelist_nemo-oce_GENERIC \
rebuild.sh \
run_nemo.sh
do

  # if a file exists with $CONFIG as suffix, then use it, otherwise take default script
  if [ -f ${TEMP_NEMO_DIR}/template_run/${file}_${CONFIG} ]; then
    if [ $file == 'run_nemo.sh' ]; then
      sed -e "s#<config>#${CONFIG}#g ; s#<case>#${CASE}#g" ${TEMP_NEMO_DIR}/template_run/${file}_${CONFIG} > run/nemo_${CONFIG}_${CASE}/${file}
      chmod +x run/nemo_${CONFIG}_${CASE}/${file}
    else
      cp -p ${TEMP_NEMO_DIR}/template_run/${file}_${CONFIG} run/nemo_${CONFIG}_${CASE}/${file}
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

# Put correct batch header :
cat ${TEMP_NEMO_DIR}/template_run/header_$HOST > tmp
cat run/nemo_${CONFIG}_${CASE}/run_nemo.sh |sed -e "s/<HEADER>//g" >> tmp
mv tmp run/nemo_${CONFIG}_${CASE}/run_nemo.sh

#mv run/nemo_${CONFIG}_${CASE}/namelist_nemo-oce_GENERIC run/nemo_${CONFIG}_${CASE}/namelist_nemo-oce_GENERIC_${CONFIG}
#mv run/nemo_${CONFIG}_${CASE}/namelist_nemo-ice_GENERIC run/nemo_${CONFIG}_${CASE}/namelist_nemo-ice_GENERIC_${CONFIG}

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

