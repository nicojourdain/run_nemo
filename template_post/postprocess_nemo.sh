#!/bin/bash
#SBATCH -C HSW24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --threads-per-core=1
#SBATCH -J post_<config>_<case>
#SBATCH -e post.e%j
#SBATCH -o post.o%j
#SBATCH --mem=32GB
#SBATCH --time=23:59:00

date

##################################################################
# USER'S CHOICES :

YEARi=1990
YEARf=2010

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

if [ ! -d ${WORKDIR}/MONTHLY ]; then
  mkdir ${WORKDIR}/MONTHLY
fi

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

echo " "
echo "========================="
date
echo "$GRID :"
echo " "

COUNT=0
NFILES=`ls -1 ${CONFIG}-${CASE}_1d_${YEAR}????_${YEAR}????_${GRID}.nc |wc -l`

if [ $NFILES -eq 0 ]; then

  echo "**************************"
  echo "* No file to postprocess *"
  echo "**************************"

else

# Check time steps consistency:
for file in ${CONFIG}-${CASE}_1d_${YEAR}????_${YEAR}????_${GRID}.nc
do
  NSTEPS=`ncdump -h $file |grep UNLIMITED |awk '{print $6}' |sed -e "s/(//g"`
  echo "${NSTEPS} time steps in $file"
  COUNT=`expr $COUNT + $NSTEPS`
done
if [ $COUNT -lt 365 ] || [ $COUNT -gt 366 ]; then
  echo "~!@#$%^&* THE TOTAL NUMBER OF TIME STEPS FOR ${GRID} AND YEAR ${YEAR} IS NOT CORRECT >>>>>> STOP !!"
  exit
fi

echo " "
# Concatenation to yearly file if needed:
if [ $NFILES -gt 1 ]; then
  rm -f tmpx_${GRID}_${YEAR}.nc
  ncrcat ${CONFIG}-${CASE}_1d_${YEAR}????_${YEAR}????_${GRID}.nc tmpx_${GRID}_${YEAR}.nc
  if [ -f tmpx_${GRID}_${YEAR}.nc ]; then
    rm -f ${CONFIG}-${CASE}_1d_${YEAR}????_${YEAR}????_${GRID}.nc
    mv tmpx_${GRID}_${YEAR}.nc ${CONFIG}-${CASE}_1d_${YEAR}0101_${YEAR}1231_${GRID}.nc
    echo "All files concatenated into ${CONFIG}-${CASE}_1d_${YEAR}0101_${YEAR}1231_${GRID}.nc"
  else
    echo "~!@#$%^&* CONCATENATION OF ${CONFIG}-${CASE}_1d_${YEAR}????_${YEAR}????_${GRID}.nc HAS FAILED >>>>> STOP !!"
    exit
  fi
else
  echo "Output file is ${CONFIG}-${CASE}_1d_${YEAR}0101_${YEAR}1231_${GRID}.nc => no need for concatenation."
fi

echo " "
# Extraction of monthly means:
if [ $COUNT == 365 ]; then
  BIS=0
else
  BIS=1
fi
DAYi=1
for MONTH in $(seq 1 12)
do
  if [ $MONTH -lt 10 ]; then
    MONTH="0${MONTH}"
  fi
  if [ $MONTH -eq 4 ] || [ $MONTH -eq 6 ] || [ $MONTH -eq 9 ] || [ $MONTH -eq 11 ]; then
    NDAYS=30
  elif [ $MONTH -eq 2 ]; then
    NDAYS=`expr 28 + $BIS`
  else
    NDAYS=31
  fi
  DAYf=`echo "${DAYi} + ${NDAYS} - 1" |bc`
  echo "$NDAYS days (${DAYi}:${DAYf}) for month ${MONTH}"
  ncks -F -d time_counter,${DAYi},${DAYf} ${CONFIG}-${CASE}_1d_${YEAR}0101_${YEAR}1231_${GRID}.nc cyp_${YEAR}_${GRID}.nc
  ncra -F -d time_counter,1,${NDAYS} cyp_${YEAR}_${GRID}.nc ../MONTHLY/${CONFIG}-${CASE}_1d_${YEAR}${MONTH}_${GRID}.nc
  rm -f cyp_${YEAR}_${GRID}.nc
  #ncra -F -d time_counter,${DAYi},${DAYf} ${CONFIG}-${CASE}_1d_${YEAR}0101_${YEAR}1231_${GRID}.nc ../MONTHLY/${CONFIG}-${CASE}_1d_${YEAR}${MONTH}_${GRID}.nc
  DAYi=`expr $DAYf + 1`
done ## MONTH

ncrcat ../MONTHLY/${CONFIG}-${CASE}_1d_${YEAR}??_${GRID}.nc ../MONTHLY/${CONFIG}-${CASE}_monthly_${YEAR}_${GRID}.nc
if [ -f ../MONTHLY/${CONFIG}-${CASE}_monthly_${YEAR}_${GRID}.nc ]; then
  rm -f ../MONTHLY/${CONFIG}-${CASE}_1d_${YEAR}??_${GRID}.nc
  echo " "
  echo "${CONFIG}-${CASE}_monthly_${YEAR}_${GRID}.nc  has been created"
else
  echo "~!@#$%^&* Concatenation of monthly file MONTHLY/${CONFIG}-${CASE}_monthly_${YEAR}_${GRID}.nc has failed >>>>> STOP !!"
  exit
fi

fi ## if [ $NFILES -eq 0 ]

done ## GRID

cd ${WORKDIR}

done ## YEAR

date
