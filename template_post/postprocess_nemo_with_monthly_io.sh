#!/bin/bash
#SBATCH -C HSW24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --threads-per-core=1
#SBATCH -J post_<config>_<case>
#SBATCH -e post.e%j
#SBATCH -o post.o%j
#SBATCH --mem=50GB
###SBATCH --time=15:29:00
#SBATCH --time=02:59:00

## FOR AMUXL12, takes ~90min per year of outputs (15h for 10 years).
## For AMUXL025, takes 4h30 to process 10 years

date

##################################################################
# USER'S CHOICES :

YEARi=2017
YEARf=2018

###################################################################

WORKDIR=`pwd`

CONFIG=$(basename $(pwd) | cut -d '_' -f2)
CASE=$(basename $(pwd) | cut -d '_' -f3-)
echo $CONFIG $CASE

for YEAR in $(seq ${YEARi} ${YEARf})
do

if [ $YEAR -lt 10 ]; then
  YEAR="0${YEAR}"
fi
if [ $YEAR -lt 100 ]; then
  YEAR="0${YEAR}"
fi
if [ $YEAR -lt 1000 ]; then
  YEAR="0${YEAR}"
fi

echo " "
echo "*********************************************************************************"
echo "*********************************************************************************"
date
echo "$CONFIG $CASE"
echo "$YEAR :"
echo " "

cd ${WORKDIR}/$YEAR

tar cvf diag_simu_${YEAR}.tar ocean.output.* namelist* app.copy 
if [ -f diag_simu_${YEAR}.tar ]; then
  rm -f ocean.output.* namelist* app.copy
else
  echo "WARNING: Something has gone wrong with creation of diag_simu_${YEAR}.tar"
  echo "  >> continue anyway..."
fi

for GRID in SBC icemod grid_T grid_U grid_V grid_W trendT trendS
do

for FREQ in 1d 5d 1m
do

echo " "
echo "========================="
date
echo "$GRID :"
echo " "

COUNT=0
NFILES=`ls -1 ${CONFIG}-${CASE}_${FREQ}_${YEAR}????_${YEAR}????_${GRID}.nc |wc -l`

if [ $NFILES -eq 0 ]; then

  echo "**************************"
  echo "* No file to postprocess *"
  echo "**************************"

else

# Check time steps consistency:
if [ $FREQ == '1d' ]; then
  for file in ${CONFIG}-${CASE}_${FREQ}_${YEAR}????_${YEAR}????_${GRID}.nc
  do
    NSTEPS=`ncdump -h $file |grep UNLIMITED |awk '{print $6}' |sed -e "s/(//g"`
    echo "${NSTEPS} time steps in $file"
    COUNT=`expr $COUNT + $NSTEPS`
  done
  if [ $COUNT -lt 365 ] || [ $COUNT -gt 366 ]; then
    echo "~!@#$%^&* THE TOTAL NUMBER OF TIME STEPS FOR ${FREQ} ${GRID} AND YEAR ${YEAR} IS NOT CORRECT >>>>>> STOP !!"
    exit
  fi
elif [ $FREQ == '1m' ]; then
  for file in ${CONFIG}-${CASE}_${FREQ}_${YEAR}????_${YEAR}????_${GRID}.nc
  do
    NSTEPS=`ncdump -h $file |grep UNLIMITED |awk '{print $6}' |sed -e "s/(//g"`
    echo "${NSTEPS} time steps in $file"
    COUNT=`expr $COUNT + $NSTEPS`
  done
  if [ $COUNT -ne 12 ]; then
    echo "~!@#$%^&* THE TOTAL NUMBER OF TIME STEPS FOR ${FREQ} ${GRID} AND YEAR ${YEAR} IS NOT CORRECT >>>>>> STOP !!"
    exit
  fi
fi

echo " "
# Concatenation to yearly file if needed:
if [ $NFILES -gt 1 ]; then
  rm -f tmpx_${FREQ}_${GRID}_${YEAR}.nc
  ncrcat ${CONFIG}-${CASE}_${FREQ}_${YEAR}????_${YEAR}????_${GRID}.nc tmpx_${FREQ}_${GRID}_${YEAR}.nc
  if [ -f tmpx_${FREQ}_${GRID}_${YEAR}.nc ]; then
    rm -f ${CONFIG}-${CASE}_${FREQ}_${YEAR}????_${YEAR}????_${GRID}.nc
    mv tmpx_${FREQ}_${GRID}_${YEAR}.nc ${CONFIG}-${CASE}_${FREQ}_${YEAR}0101_${YEAR}1231_${GRID}.nc
    echo "All files concatenated into ${CONFIG}-${CASE}_${FREQ}_${YEAR}0101_${YEAR}1231_${GRID}.nc"
  else
    echo "~!@#$%^&* CONCATENATION OF ${CONFIG}-${CASE}_${FREQ}_${YEAR}????_${YEAR}????_${GRID}.nc HAS FAILED >>>>> STOP !!"
    exit
  fi
else
  echo "Output file is ${CONFIG}-${CASE}_${FREQ}_${YEAR}0101_${YEAR}1231_${GRID}.nc => no need for concatenation."
fi

fi ## if [ $NFILES -eq 0 ]

done ## $FREQ

done ## GRID

cd ${WORKDIR}

done ## YEAR

date
