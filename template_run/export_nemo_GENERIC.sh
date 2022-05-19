<HEADER>

set -x
ulimit -s unlimited

#- Netcdf library for small fortran scripts (not for NEMO)
export NC_INC="-I`nc-config --includedir`"
export NC_LIB=`nc-config --flibs`

WORKDIR='WWWW'
NRUN=NNNN
YEAR='YYYY'
STOCKDIR='SSSS'
CONFIG='CCCC'
CONFPAR='cccc'
CASE='OOOO'
NZOOM=ZZZZ
NITEND=UUUU

YEARp1=`expr $YEAR + 1`
if [ $YEARp1 -lt 1000 ]; then
  YEARp1="0$YEARp1"
fi

echo "****************************************"
echo "*   Netcdf compression then export      "
echo "****************************************"

date

###################################################################
##-- rebuild NEMO outputs (nests + parent grid)
##   if xios's iodef.xml is in "multiple_file" mode :

cd $WORKDIR/OUTPUT_${NRUN}

MULTIPLE_FILE=`grep multiple_file ../file_def_nemo-oce.xml  | wc -l`

if [ $MULTIPLE_FILE -gt 0 ]; then

    echo "############################################################################"
    echo "## xios is in multiple_file mode ==> rebuilding outputs into single files ##"

    for file in `ls -1 *${CONFPAR}-${CASE}_*_0000.nc`
    do

      prefix_file=`echo $file | sed -e "s/_0000\.nc//g"`
      sed -e "s/pppppp/${prefix_file}/g" ../rebuild_outputs_GENERIC.f90 > rebuild_outputs.f90
      ifort -c $NC_INC rebuild_outputs.f90
      ifort -o rebuild_outputs rebuild_outputs.o $NC_LIB
      ./rebuild_outputs
      rm -f rebuild_outputs rebuild_outputs.o rebuild_outputs.f90

      file_out=`echo $file | sed -e "s/_0000//g"`
      if [ -f ${file_out} ]; then
        echo "${file_out}   [oK]"
        rm -f ${prefix_file}_[0-9][0-9][0-9][0-9].nc
      else
        echo "~!@#%^&* ERROR : ${file_out} HAS NOT BEEN CREATED !!!"
        echo "                 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> STOP !!"
        exit
      fi

    done

else

  echo "## xios is in single file mode ==> no need to rebuild"

fi

##########################################################
##-- compress output files:
##   ( read list of files from this run only, not to read 
##     outputs from new on going simulation )

echo " "
ls -1 *_[1-5][d-y]_*.nc
echo " "

rm -f tmp.nc

#NFILES=`ls -1 *_[1-5][d-m]_*.nc | wc -l | cut -d ' ' -f1`
#ifile=1
#while [ $ifile -le $NFILES ];
#do
#   file=`ls -1 *_[1-5][d-m]_*.nc | head -${ifile} | tail -1`
#   nccopy -d 9 $file tmp.nc
#   if [ -f tmp.nc ]; then
#      mv -f tmp.nc $file
#   else
#      echo "!@#$% PROBLEM WITH COMPRESSION OF FILE $file"
#      echo " >>>>>>>>>> STOP !"
#      exit
#   fi 
#   ifile=`expr $ifile + 1`
#done

cd $WORKDIR

##########################################################
##-- Prepare directories

if [ ! -d ${STOCKDIR}/output ]; then
  mkdir ${STOCKDIR}/output
fi

if [ ! -d ${STOCKDIR}/restart ]; then
  mkdir ${STOCKDIR}/restart
fi

if [ ! -d ${STOCKDIR}/output/nemo_${CONFIG}_${CASE} ]; then
  mkdir ${STOCKDIR}/output/nemo_${CONFIG}_${CASE}
fi

if [ ! -d ${STOCKDIR}/output/nemo_${CONFIG}_${CASE}/${YEAR} ]; then
  mkdir ${STOCKDIR}/output/nemo_${CONFIG}_${CASE}/${YEAR}
fi

if [ ! -d ${STOCKDIR}/output/nemo_${CONFIG}_${CASE}/${YEARp1} ]; then
  mkdir ${STOCKDIR}/output/nemo_${CONFIG}_${CASE}/${YEARp1}
fi

if [ ! -d ${STOCKDIR}/restart/nemo_${CONFIG}_${CASE} ]; then
  mkdir ${STOCKDIR}/restart/nemo_${CONFIG}_${CASE}
fi

#############################################
#- Output files

echo " "
echo "Moving output files to ${STOCKDIR}..."

mv -f OUTPUT_${NRUN}/*  ${STOCKDIR}/output/nemo_${CONFIG}_${CASE}/${YEAR} 

#############################################
#- Restart files

# NB: starting to copy last restart date (for next re-submission)

echo " "
echo "Moving restart files to ${STOCKDIR}..."

## NITEND with zeros in front:
ENDLONG=`echo "${NITEND} + 100000000" |bc | sed -e "s/1//"`

NRST=`ls -1 ${WORKDIR}/restart_????????.nc | wc -l`
for iRST in $(seq 2 $NRST) ## NB: the last restart file is kept for the next run
do
  DTMP_RST=`ls -1 restart_????????.nc | tail -${iRST} | head -1`
  echo ${DTMP_RST}
  NIT_RST=`basename ${DTMP_RST} | cut -d '_' -f2 | cut -d '.' -f1`
  mv restart_${NIT_RST}.nc ${STOCKDIR}/restart/nemo_${CONFIG}_${CASE}/.
done

NRST=`ls -1 ${WORKDIR}/restart_ice_????????.nc | wc -l`
for iRST in $(seq 2 $NRST) ## NB: the last restart file is kept for the next run
do
  DTMP_RST=`ls -1 restart_ice_????????.nc | tail -${iRST} | head -1`
  echo ${DTMP_RST}
  NIT_RST=`basename ${DTMP_RST} | cut -d '_' -f3 | cut -d '.' -f1`
  mv restart_ice_${NIT_RST}.nc ${STOCKDIR}/restart/nemo_${CONFIG}_${CASE}/.
done

for iZOOM in $(seq 1 $NZOOM)
do
  NRST=`ls -1 ${WORKDIR}/${iZOOM}_restart_????????.nc | wc -l`
  for iRST in $(seq 2 $NRST) ## NB: the last restart file is kept for the next run
  do
    DTMP_RST=`ls -1 ${iZOOM}_restart_????????.nc | tail -${iRST} | head -1`
    echo ${DTMP_RST}
    NIT_RST=`basename ${DTMP_RST} | cut -d '_' -f3 | cut -d '.' -f1`
    mv ${iZOOM}_restart_${NIT_RST}.nc ${STOCKDIR}/restart/nemo_${CONFIG}_${CASE}/.
  done
  NRST=`ls -1 ${WORKDIR}/${iZOOM}_restart_ice_????????.nc | wc -l`
  for iRST in $(seq 2 $NRST) ## NB: the last restart file is kept for the next run
  do
    DTMP_RST=`ls -1 ${iZOOM}_restart_ice_????????.nc | tail -${iRST} | head -1`
    echo ${DTMP_RST}
    NIT_RST=`basename ${DTMP_RST} | cut -d '_' -f4 | cut -d '.' -f1`
    mv ${iZOOM}_restart_ice_${NIT_RST}.nc ${STOCKDIR}/restart/nemo_${CONFIG}_${CASE}/.
  done
done

echo "Transfer completed."
echo " "
date
