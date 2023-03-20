<HEADER>

date

set -x
ulimit -s unlimited

#=================================================================================
#=================================================================================
# 0- User's choices
#=================================================================================
#=================================================================================

MACHINE='<host>' ## Either 'occigen' or 'irene'

CONFIG='<config>'  ## FULL CONFIG NAME (e.g. "AMUXL025.L75", "trop075", "trop075_nest025")
                       ## NB: THIS NAME SHOULD NOT START WITH A NUMBER

CONFPAR=$CONFIG #- IF NO NEST SHOULD BE EQUAL TO $CONFIG
                #  IF NESTS, SHOULD BE THE ABSOLUTE PARENT CONFIG NAME
                #  (e.g. CONFPAR="trop075" when CONFIG="trop075_nest025")

CONFEXE='WED025'   # only for nemo.exe

CASE='<case>' # should not be too long [>15 char.] otherwise, NEMO file names are affected

YEAR0=1979      #- initial year of the long experiment  [ 4 digits ]
MONTH0=03       #- initial month of the long experiment [ 2 digits ]
DAY0=01         #- initial day of the long experiment   [ 2 digits ]

YEAR_MAX=2018   #- stop after $YEAR_MAX is completed

NRUN_MAX=2      #- stop after $NRUN_MAX re-submissions

NDAYS=1         #- Split year by slices of $NDAYS days (possibly less if BYMONTH=1)

BYMONTH=0       # =0 to stick to a duration of NDAYS (or cut to either 365 or 366 days if > 1year).
                # =1 to cut runs to entire months (e.g. cut to 31 days for January if NDAYS=32).

OUTPUT_FREQ='1d' # = '1d' for daily-mean outputs
                 # = '5d' for 5-day-mean outputs (only for BYMONTH=0)
                 # = '1m' for monthly-mean outputs (only for BYMONTH=1)
                 # = '1d1m' for both daily-mean and monthly-mean outputs (only for BYMONTH=1)
                 # = '5d1y' for both 5-day-mean and yearly-mean outputs (only for BYMONTH=0 and full years) -- check for leap years
                 # = '1m1y' for both monthly-mean and yearly-mean outputs (only for BYMONTH=1 and full years)
                 # = '5d1m1y' for both 5-day-mean and yearly-mean outputs (only for BYMONTH=1 and full years) -- check for leap years

WORKDIR=`pwd`

STOCKDIR="$GEN6035_CCCSCRATCHDIR"  #- restart, output directory

INPUTDIR="${GEN6035_CCCWORKDIR}/input/nemo_${CONFIG}"  #- input directory

GLOBAL_SIM='OPM026'  # ORCA simulation used for BDY, runoff and sss restoring
                     # (also used to look for possible istate file names)
                     # Define GLOBAL_SIM='' if not used.
BDYDIR="${INPUTDIR}/BDY_${GLOBAL_SIM}"  #- input directory for BDYs
SSSDIR="${INPUTDIR}/SSS_${GLOBAL_SIM}"  #- input directory for SSS relaxation (if any)
RNFDIR="${INPUTDIR}/RNF_${GLOBAL_SIM}"  #- input directory for runoff relaxation (if any)

TOPO="BedMachineAntarctica-2020-10-08"  # to use domain_cfg_${CONFPAR}_${TOPO}.nc
                                        # define TOPO="" to use domain_cfg_${CONFPAR}.nc

FORCINGdir="${GEN6035_ALL_CCCWORKDIR}/MODEL_INPUTS/NEMO/FORCINGS_SETS/JRA55" # Atmospheric forcing

#- Netcdf library for small fortran scripts (not for NEMO)
export NC_INC="-I`nc-config --includedir`"
export NC_LIB=`nc-config --flibs`

NEMOdir="${GEN6035_CCCWORKDIR}/models/nemo_4.2.0" # NEMO model directory
XIOSdir="${GEN6035_CCCWORKDIR}/models/xios_trunk" # XIOS directory

NZOOM=0  # nb of agrif nests (0 if no agrif nest)

NB_TASK_PER_NODE=128    # 128 (TGCC), 28 (OCCIGEN)
NB_TASK_XIOS_PER_NODE=4 # Number of core used per xios on each node (should typically be in the 1-3 range).

#=================================================================================
#=================================================================================
# 1- Initialization
#=================================================================================
#=================================================================================

PWDDIR=`pwd`
if [ ! `basename $PWDDIR` == nemo_${CONFIG}_${CASE} ]; then
 echo '~!@#%^&* CHECK CONFIG and CASE in run_nemo.sh >>>>>>>>>> Stop !!'
 date
 exit
fi

export NB_NODES=`echo "$SLURM_NTASKS / $NB_TASK_PER_NODE" |bc`
export NB_TASK_NEMO_PER_NODE=$(( NB_TASK_PER_NODE - NB_TASK_XIOS_PER_NODE ))
export NB_TASK_XIOS=`echo "$NB_NODES * $NB_TASK_XIOS_PER_NODE" |bc`
export NB_TASK_NEMO=$(( SLURM_NTASKS - NB_TASK_XIOS ))

## { unset initiaux 
#unset    OMPI_MCA_ess
##
#unset    OMPI_MCA_pml
#unset    OMPI_MCA_mtl
#unset    OMPI_MCA_mtl_mxm_np 
#unset    OMPI_MCA_pubsub  
## }

############################################################
##-- create links to executables :

rm -f nemo.exe
if [ $NZOOM -gt 0 ]; then
  ln -s ${NEMOdir}/cfgs/${CONFEXE}_agrif/BLD/bin/nemo.exe
else
  ln -s ${NEMOdir}/cfgs/${CONFEXE}/BLD/bin/nemo.exe
fi

rm -f xios_server.exe
ln -s ${XIOSdir}/bin/xios_server.exe

rm -f file_def_nemo-oce.xml file_def_nemo-ice.xml
if [ $BYMONTH -eq 1 ]; then
  if [ $OUTPUT_FREQ == '1d1m' ] || [ $OUTPUT_FREQ == '1m' ] || [ $OUTPUT_FREQ == '1d' ]; then
    ln -s -v file_def_nemo-oce_${OUTPUT_FREQ}.xml file_def_nemo-oce.xml
    ln -s -v file_def_nemo-ice_${OUTPUT_FREQ}.xml file_def_nemo-ice.xml
  else
    echo '~!@#%^&* wrong value of OUTPUT_FREQ >>>>>>>>>> Stop !!'
    exit
  fi
else
  if [ $OUTPUT_FREQ == '1d' ] || [ $OUTPUT_FREQ == '5d' ] || [ $OUTPUT_FREQ == '5d1y' ] || [ $OUTPUT_FREQ == '5d1m1y' ]; then
    ln -s -v file_def_nemo-oce_${OUTPUT_FREQ}.xml file_def_nemo-oce.xml
    ln -s -v file_def_nemo-ice_${OUTPUT_FREQ}.xml file_def_nemo-ice.xml
  else
    echo '~!@#%^&* wrong value of OUTPUT_FREQ >>>>>>>>>> Stop !!'
    exit
  fi
fi

##############################################
##-- define current year and nb of days

if [ $BYMONTH -eq 1 ] && [ $NDAYS -lt 28 ]; then
  echo "~!@#% You need NDAYS > 27 if BYMONTH = 1"
  echo "             >>>>> STOP"
  date
  exit
fi

if [ -f prod_nemo.db ]; then
read NRUN YEAR MONTH DAY NITENDM1 NITENDM1ZOOM << EOF
`tail -1 prod_nemo.db`
EOF
else
echo "1 ${YEAR0} ${MONTH0} ${DAY0} 0" > prod_nemo.db
YEAR=${YEAR0}
MONTH=${MONTH0}
DAY=${DAY0}
NRUN=1
NITENDM1=0  ## last time step of previous run
## add last time step of previous runs on children domains 
## at the end of the line in prod_nemo.db :
for iZOOM in $(seq 1 ${NZOOM}); do
  sed -e "s/$/ 0/g" prod_nemo.db > tmp
  mv tmp prod_nemo.db
done
fi

if [ $YEAR -gt ${YEAR_MAX} ]; then
  echo " "
  echo "Year greater then YEAR_MAX >>>>>>>> stop !!"
  date
  exit
fi

if [ $NRUN -gt ${NRUN_MAX} ]; then
  echo " "
  echo "Nb of re-submission greater then NRUN_MAX >>>>>>>> stop !!"
  date
  exit
fi

#####
# adjust nb of days to finish at the end of the year
ISLEAP=`grep nn_leapy namelist_nemo-oce_GENERIC | awk '{print $3}'`
if [ $BYMONTH -eq 0 ]; then
if [ ! -f calculate_end_date ]; then
  ifort -o calculate_end_date calculate_end_date.f90
fi
read YEARf MONTHf DAYf NDAYScorr << EOF
`./calculate_end_date $YEAR $MONTH $DAY $NDAYS $ISLEAP`
EOF
if [ $NDAYScorr -ne $NDAYS ]; then
echo "Adjusting run length to finish at the end of current year"
NDAYS=$NDAYScorr
read YEARf MONTHf DAYf NDAYScorr << EOF
`./calculate_end_date $YEAR $MONTH $DAY $NDAYS $ISLEAP`
EOF
fi
elif [ $BYMONTH -eq 1 ]; then
if [ ! -f calculate_end_date_month ]; then
  ifort -o calculate_end_date_month calculate_end_date_month.f90
fi
read YEARf MONTHf DAYf NDAYScorr << EOF
`./calculate_end_date_month $YEAR $MONTH $DAY $NDAYS $ISLEAP`
EOF
if [ $NDAYScorr -ne $NDAYS ]; then
echo "Adjusting run length to finish at the end of a month"
NDAYS=$NDAYScorr
read YEARf MONTHf DAYf NDAYScorr << EOF
`./calculate_end_date_month $YEAR $MONTH $DAY $NDAYS $ISLEAP`
EOF
fi
else
  echo "~!@@#$%^& WRONG VALUE FOR 'BYMONTH'  >>>> stop"
  date
  exit
fi

if [ $NDAYS -gt 0 ]; then
  echo " -> Run duration = $NDAYS days"
  echo " "
else
  echo " -> Run duration < 1   --> STOP !@#%^&"
  echo " (check BYMONTH option)"
  date
  exit
fi

##-- calculate corresponding number of time steps for NEMO:
RN_DT=`grep "rn_Dt " namelist_nemo-oce_GENERIC |grep 'dynamics' |grep 'step' |cut -d '=' -f2 | cut -d '!' -f1 | sed -e "s/ //g"`
NIT000=`echo "$NITENDM1 + 1" | bc`
NITEND=`echo "$NITENDM1 + ${NDAYS} * 86400 / ${RN_DT}" | bc`

echo "****************************************************"
echo "*          NEMO SIMULATION                          "
echo "*   config  $CONFIG                                 "
echo "*   case    $CASE                                   "
echo "*   from    ${DAY}/${MONTH}/${YEAR}                 "
echo "*   to      ${DAYf}/${MONTHf}/${YEARf}              "
echo "*   i.e. step $NIT000 to $NITEND (for mother grid)  "
echo "*                                                   "
echo "*   total number of tasks >>>>> ${SLURM_NTASKS}     "
echo "*   number of xios tasks  >>>>> ${NB_TASK_XIOS}     "
echo "*                                                   "
echo "****************************************************"
echo " "
date
echo " "

#####################################################################
##-- create executable and rebuild namelist to rebuild restart files

###############################################################
##-- edit NEMO's namelist

echo "Editing namelist..."

rm -f namelist_ref namelist_cfg

##- Start from restart file or at rest from initial T,S conditions:
if [ $NRUN -gt 1 ] || [ $CONFIG == "trop075" ]; then
  sed -e "s#<RESTNEM>#.true.#g"  namelist_nemo-oce_GENERIC > namelist_ref
  RST=1
else
  sed -e "s#<RESTNEM>#.false.#g" namelist_nemo-oce_GENERIC > namelist_ref
  RST=0
fi

##- Specific treatment for TROP075's restart/initial state:
if [ $NRUN -eq 1 ] && [ $CONFIG == "trop075" ]; then
  sed -e "s#<AAAA>#  0 #g" namelist_ref > tmp
  mv -f tmp namelist_ref
else
  sed -e "s#<AAAA>#  2 #g"  namelist_ref > tmp
  mv -f tmp namelist_ref
fi

##- Ice Sheet Coupling:
IS_ISCPL=`grep ln_isfcpl namelist_nemo-oce_GENERIC | grep -v ln_isfcpl_cons | awk '{print $3}' |sed -e "s/\.//g"`
if [ $IS_ISCPL == 'true' ]; then
  echo "WARNING : enabling the ice shelf geometry to move (ln_isfcpl=true) !!!"
fi

echo "s#<CCCC>#${CONFIG}#g ; s#<OOOO>#${CASE}#g ; s#<IIII>#${YEAR0}${MONTH0}${DAY0}#g ; s#<NIT000>#${NIT000}#g ; s#<NITEND>#${NITEND}#g"
sed -e "s#<CCCC>#${CONFIG}#g ; s#<OOOO>#${CASE}#g ; s#<IIII>#${YEAR0}${MONTH0}${DAY0}#g ; s#<NIT000>#${NIT000}#g ; s#<NITEND>#${NITEND}#g" namelist_ref > tmp
mv -f tmp namelist_ref

ln -s namelist_ref namelist_cfg

for iZOOM in $(seq 1 $NZOOM); do
  echo "Editing ${iZOOM}_namelist..."
  rm -f ${iZOOM}_namelist_ref ${iZOOM}_namelist_cfg
  if [ $NRUN -gt 1 ] || [ $CONFIG == "trop075" ]; then
    sed -e "s#<RESTNEM>#.true.#g"  ${iZOOM}_namelist_nemo-oce_GENERIC > ${iZOOM}_namelist_ref
    RST=1
  else
    sed -e "s#<RESTNEM>#.false.#g" ${iZOOM}_namelist_nemo-oce_GENERIC > ${iZOOM}_namelist_ref
    RST=0
  fi
  ##- Specific treatment for TROP075's restart/initial state:
  if [ $NRUN -eq 1 ] && [ $CONFIG == "trop075" ]; then
    sed -e "s#<AAAA>#  0 #g" ${iZOOM}_namelist_ref > tmp
    mv -f tmp ${iZOOM}_namelist_ref
  else
    sed -e "s#<AAAA>#  2 #g"  ${iZOOM}_namelist_ref > tmp
    mv -f tmp ${iZOOM}_namelist_ref
  fi
  ##- calculate initial and last time step for the child domains :
  ##-- calculate corresponding number of time steps for NEMO:
  RN_DT_ZOOM=`grep "rn_Dt " ${iZOOM}_namelist_nemo-oce_GENERIC |cut -d '=' -f2 | cut -d '!' -f1 | sed -e "s/ //g"`
  NIT000_ZOOM=`echo "( ${NITENDM1} * ${RN_DT} / ${RN_DT_ZOOM} ) + 1" | bc`
  NITEND_ZOOM=`echo "( ${NITENDM1} * ${RN_DT} / ${RN_DT_ZOOM} ) + ${NDAYS} * 86400 / ${RN_DT_ZOOM}" | bc`
  ##--
  sed -e "s#<CCCC>#${CONFIG}#g ; s#<OOOO>#${CASE}#g ; s#<IIII>#${YEAR0}0101#g ; s#<NIT000>#${NIT000_ZOOM}#g ; s#<NITEND>#${NITEND_ZOOM}#g" ${iZOOM}_namelist_ref > tmp
  mv -f tmp ${iZOOM}_namelist_ref
  ln -s ${iZOOM}_namelist_ref ${iZOOM}_namelist_cfg
done

rm -f namelist_ice_ref namelist_ice_cfg
cp -p namelist_nemo-ice_GENERIC namelist_ice_ref
ln -s namelist_ice_ref namelist_ice_cfg

#############################################################
###-- prepare script that will be used to export outputs :

STOCKDIR2=`echo $STOCKDIR |sed -e "s/\//\\\\\\\\\//g"`
WORKDIR2=`echo $WORKDIR  |sed -e "s/\//\\\\\\\\\//g"`

sed -e "s/CCCC/${CONFIG}/g ; s/cccc/${CONFPAR}/g ; s/OOOO/${CASE}/g ; s/SSSS/${STOCKDIR2}/g ; s/WWWW/${WORKDIR2}/g ; s/YYYY/${YEAR}/g ; s/NNNN/${NRUN}/g ; s/ZZZZ/${NZOOM}/g ; s/UUUU/${NITEND}/g" export_nemo_GENERIC.sh > export_nemo_${NRUN}.sh

chmod +x export_nemo_${NRUN}.sh

#=======================================================================================
#=======================================================================================
# 2- Manage links to input files
#=======================================================================================
#=======================================================================================

echo " "
date
echo " "
echo " Linking input files from ${INPUTDIR}"

DATE0=`grep nn_date0 namelist_cfg | head -1 | awk '{print $3}'`
Y0=`echo $DATE0 | cut -c 1-4`
M0=`echo $DATE0 | cut -c 5-6`

YEARm1=`expr $YEAR - 1`
if [ $YEARm1 -lt 1000 ]; then
  YEARm1="0$YEARm1"
fi
# because no data before:
if [ $YEAR -eq $Y0 ]; then
  YEARm1=$Y0
fi
YEARp1=`expr $YEAR + 1`
if [ $YEARp1 -lt 1000 ]; then
  YEARp1="0$YEARp1"
fi

#####################################################################
##-- import files that are not time dependent if not already there :

# Function that tells if namelist variable (1st argument) is in 
# a given namelist section (2nd arg.) for a given namelist (3rd arg.):
get_value_in_namelist() {
  sed -n "/^${2}/, /^\//p" $3 |grep $1 | cut -d '=' -f2 | cut -d '!' -f1 | sed -e "s/ //g"
}

# Loop over parent config (iZOOM=0) and AGRIF child domains (iZOOM>0) :
for iZOOM in $(seq 0 ${NZOOM}); do

  if [ $iZOOM -eq 0 ]; then
    PREFIX=''
  else
    PREFIX="${iZOOM}_"
  fi

  ##-- All geographical and grid characteristics:
  rm -f ${PREFIX}domain_cfg.nc
  if [ -z $TOPO ]; then
    ln -s -v ${INPUTDIR}/${PREFIX}domain_cfg_${CONFPAR}.nc ${PREFIX}domain_cfg.nc
  else
    ln -s -v ${INPUTDIR}/${PREFIX}domain_cfg_${CONFPAR}_${TOPO}.nc ${PREFIX}domain_cfg.nc
  fi
  
  ##-- Chlorophyll Concentration (for vertical penetration of solar flux):
  IS_TRAQSR=`get_value_in_namelist 'ln_traqsr' '&namsbc' ${PREFIX}namelist_cfg`
  if [ $IS_TRAQSR == ".true." ]; then
    rm -f ${PREFIX}chlorophyll.nc
    if [ -f ${INPUTDIR}/${PREFIX}chlorophyll_${CONFPAR}.nc ]; then
      # use interpolated file :
      ln -s -v ${INPUTDIR}/${PREFIX}chlorophyll_${CONFPAR}.nc ${PREFIX}chlorophyll.nc
    else
      # use weights :
      WGT_QSR_FILE=`grep sn_chl ${PREFIX}namelist_cfg | grep '\.' | cut -d ',' -f7 |sed -e "s/'//g"`
      if [ -f ${INPUTDIR}/${PREFIX}${WGT_QSR_FILE} ]; then
        ln -s -v ${INPUTDIR}/${PREFIX}chlorophyll.nc
        ln -s -v ${INPUTDIR}/${PREFIX}${WGT_QSR_FILE}
      fi
    fi
  fi
  
  ##-- Geothermal Heat Flux :
  IS_GEO1=`get_value_in_namelist 'ln_trabbc' '&nambbc' ${PREFIX}namelist_cfg`
  IS_GEO2=`get_value_in_namelist 'nn_geoflx' '&nambbc' ${PREFIX}namelist_cfg`
  if [ $IS_GEO1 == ".true." ] && [ $IS_GEO2 == "2" ]; then
    rm -f ${PREFIX}geothermal_heating.nc
    if [ -f ${INPUTDIR}/${PREFIX}geothermal_heating_${CONFPAR}.nc ]; then
      # use interpolated file :
      ln -s -v ${INPUTDIR}/${PREFIX}geothermal_heating_${CONFPAR}.nc ${PREFIX}geothermal_heating.nc
    else
      # use weights :
      WGT_GEO_FILE=`grep sn_qgh ${PREFIX}namelist_cfg | grep '\.' | cut -d ',' -f7 |sed -e "s/'//g"`
      if [ -f ${INPUTDIR}/${PREFIX}${WGT_GEO_FILE} ]; then
        ln -s -v ${INPUTDIR}/${PREFIX}geothermal_heating.nc
        ln -s -v ${INPUTDIR}/${PREFIX}${WGT_GEO_FILE}
      fi
    fi
  fi
  
  ##-- Liquid (and solid) runoff :
  IS_RNF=`get_value_in_namelist 'ln_rnf' '&namsbc' ${PREFIX}namelist_cfg`
  if [ $IS_RNF == ".true." ]; then
    IS_CLIM=`grep " sn_rnf " ${PREFIX}namelist_cfg | cut -d ',' -f5 | sed -e "s/ //g"`
    if [ $IS_CLIM == ".true." ]; then
      rm -f ${PREFIX}runoff.nc ${PREFIX}runoff_y????.nc
      ln -s -v ${INPUTDIR}/${PREFIX}runoff_${CONFPAR}.nc ${PREFIX}runoff.nc
    else
      for AN in $YEARm1 $YEAR $YEARp1; do
        rm -f ${PREFIX}runoff.nc ${PREFIX}runoff_y${AN}.nc
        if [ -f ${RNFDIR}/${PREFIX}runoff_y${AN}_${CONFPAR}.nc ]; then 
          ln -s -v ${RNFDIR}/${PREFIX}runoff_y${AN}_${CONFPAR}.nc ${PREFIX}runoff_y${AN}.nc
        fi
      done
    fi
  fi
  
  ##-- 3D T,S restoring :
  IS_TRADMP=`get_value_in_namelist 'ln_tradmp' '&namtra_dmp' ${PREFIX}namelist_cfg`
  if [ $IS_TRADMP == ".true." ]; then
    rm -f ${PREFIX}resto.nc
    ln -s -v ${INPUTDIR}/${PREFIX}resto_${CONFPAR}.nc ${PREFIX}resto.nc
  fi
  
  ##-- Internal Wave Mixing :
  IS_IWM=`get_value_in_namelist 'ln_mevar' '&namzdf_iwm' ${PREFIX}namelist_cfg`
  if [ $IS_IWM == ".true." ]; then
   rm -f ${PREFIX}iwm.nc
   ln -s -v ${INPUTDIR}/${PREFIX}iwm_${CONFIG}.nc ${PREFIX}iwm.nc
  fi
  
  ##-- Localy-boosted Bottom Friction :
  IS_BFR_BOOST=`get_value_in_namelist 'ln_boost' '&namdrg_bot' ${PREFIX}namelist_cfg`
  if [ $IS_BFR_BOOST == ".true." ]; then
    rm -f ${PREFIX}bfr_coef.nc
    ln -s -v ${INPUTDIR}/${PREFIX}bfr_${CONFIG}.nc ${PREFIX}bfr_coef.nc
  fi
  
  ##-- Localy-boosted Top Friction :
  IS_TFR_BOOST=`get_value_in_namelist 'ln_boost' '&namdrg_top' ${PREFIX}namelist_cfg`
  if [ $IS_TFR_BOOST == ".true." ]; then
    rm -f ${PREFIX}tfr_coef.nc
    ln -s -v ${INPUTDIR}/${PREFIX}tfr_${CONFIG}.nc ${PREFIX}tfr_coef.nc
  fi
  
  ##-- Localy-varying lateral friction (DRAKKAR only) :
  IS_DRAK=`grep '&namlbc_drk' ${PREFIX}namelist_cfg |wc -l`
  if [ $IS_DRAK -ne 0 ]; then
    LS_SHLAT2D=`get_value_in_namelist 'ln_shlat2d' '&namlbc_drk' ${PREFIX}namelist_cfg`
    if [ $LS_SHLAT2D == ".true." ]; then
      rm -f ${PREFIX}shlat2d.nc
      ln -s -v ${INPUTDIR}/${PREFIX}shlat_${CONFPAR}.nc ${PREFIX}shlat2d.nc
    fi
  fi

  ##########
  ##-- import Boundary Conditions (BDYs) 

  IS_BDY=`get_value_in_namelist 'ln_bdy' '&nambdy' ${PREFIX}namelist_cfg` 
  
  if [ $IS_BDY == ".true." ]; then

    rm -f coordinates_bdy.nc
    ln -s -v ${INPUTDIR}/coordinates_bdy_${CONFPAR}.nc coordinates_bdy.nc
    
    for AN in $YEARm1 $YEAR $YEARp1; do
      for BDYNAM in bdyT_tra bdyU_u2d bdyU_u3d bdyV_u2d bdyV_u3d bdyT_ice bdyT_ssh; do
        rm -f ${PREFIX}${BDYNAM}_y${AN}.nc
        if [ -f ${BDYDIR}/${PREFIX}${BDYNAM}_y${AN}_${CONFPAR}.nc ]; then
          ln -s -v ${BDYDIR}/${PREFIX}${BDYNAM}_y${AN}_${CONFPAR}.nc ${PREFIX}${BDYNAM}_y${AN}.nc
        fi
      done
    done
    
    ## TIDES :
    for HARM in `grep sn_tide_cnames ${PREFIX}namelist_cfg | awk '{print $3}' |sed -e "s/'//g"`; do
      ln -s -v ${INPUTDIR}/BDY_TIDES/${PREFIX}bdytide_${CONFPAR}_${HARM}_grid_T.nc ${PREFIX}bdytide_${HARM}_grid_T.nc
      ln -s -v ${INPUTDIR}/BDY_TIDES/${PREFIX}bdytide_${CONFPAR}_${HARM}_grid_U.nc ${PREFIX}bdytide_${HARM}_grid_U.nc
      ln -s -v ${INPUTDIR}/BDY_TIDES/${PREFIX}bdytide_${CONFPAR}_${HARM}_grid_V.nc ${PREFIX}bdytide_${HARM}_grid_V.nc
    done

  fi # IS_BDY
  
  ################################
  ##-- atmospheric forcing
  
  FORDTA=`basename $FORCINGdir`

  rm -f ${PREFIX}w_bilin.nc ${PREFIX}w_bicubic.nc
  ln -s -v ${INPUTDIR}/${PREFIX}weights_bilin_${FORDTA}_${CONFIG}.nc  ${PREFIX}w_bilin.nc
  ln -s -v ${INPUTDIR}/${PREFIX}weights_bicub_${FORDTA}_${CONFIG}.nc  ${PREFIX}w_bicubic.nc
  
  # The following only work if in thefollowing order: namsbc_clio, namsbc_core, namsbc_mfs
  IS_BLK=`get_value_in_namelist 'ln_blk' '&namsbc' ${PREFIX}namelist_cfg`
  if [ $IS_BLK == ".true." ]; then
    LINE_BLK=`grep -n namsbc_blk ${PREFIX}namelist_cfg | tail -1 | cut -d ':' -f1`
    for NAMAT in sn_wndi sn_wndj sn_qsr sn_qlw sn_tair sn_humi sn_prec sn_snow sn_slp; do
      ATM_FILE=`awk "/${NAMAT}/ && NR >= ${LINE_BLK}" ${PREFIX}namelist_cfg | cut -d "'" -f2 | grep -v \"${NAMAT}\" | head -1`
      IS_CLIM=`awk "/${NAMAT}/ && NR >= ${LINE_BLK}" ${PREFIX}namelist_cfg | grep -v \"${NAMAT}\" | cut -d ',' -f5 | sed -e "s/ //g" | head -1`
      if [ $IS_CLIM == ".true." ]; then
        rm -f ${PREFIX}${ATM_FILE}.nc
        if [ -f ${FORCINGdir}/${PREFIX}${ATM_FILE}.nc ]; then
          ln -s -v ${FORCINGdir}/${PREFIX}${ATM_FILE}.nc
        fi
      else
        for AN in  $YEARm1 $YEAR $YEARp1; do
          rm -f ${PREFIX}${ATM_FILE}_y${AN}.nc
          if [ -f ${FORCINGdir}/${PREFIX}${ATM_FILE}_y${AN}.nc ]; then
            ln -s -v ${FORCINGdir}/${PREFIX}${ATM_FILE}_y${AN}.nc ${PREFIX}${ATM_FILE}_y${AN}.nc
          fi
        done
      fi
    done
  fi

  ##########################
  ##-- SSS restoring if any :
  
  SSSREL=`get_value_in_namelist 'nn_sssr' '&namsbc_ssr' ${PREFIX}namelist_cfg`
  
  if [ $SSSREL -gt 0 ]; then
    echo "Run with SSS relaxation : nn_sssr = $SSSREL"
    SSS_FILE=`grep sn_sss ${PREFIX}namelist_cfg | cut -d "'" -f2 | head -1`
    IS_CLIM=`grep sn_sss ${PREFIX}namelist_cfg | cut -d ',' -f5 | sed -e "s/ //g" | head -1`
    if [ $IS_CLIM == ".true." ]; then
      rm -f ${PREFIX}${SSS_FILE}.nc
      if [ -f ${INPUTDIRDIR}/${PREFIX}sss.nc ]; then
        ln -s -v ${INPUTDIRDIR}/${PREFIX}sss.nc ${PREFIX}${SSS_FILE}.nc
      fi
    else
      for AN in  $YEARm1 $YEAR $YEARp1 ; do
        rm -f ${PREFIX}${SSS_FILE}_y${AN}.nc
        if [ -f ${SSSDIR}/${PREFIX}sss_y${AN}_${CONFIG}.nc ]; then
          ln -s -v ${SSSDIR}/${PREFIX}sss_y${AN}_${CONFIG}.nc ${PREFIX}${SSS_FILE}_y${AN}.nc
        fi
      done
    fi
  fi
  
  ################################
  ##-- Initial state or Restart
  
  rm -f ${PREFIX}restart.nc ${PREFIX}restart_ice.nc
  rm -f ${PREFIX}istate_TS_y????m??.nc ${PREFIX}istate_sea_ice_y????m??.nc ${PREFIX}istate_TS.nc ${PREFIX}istate_sea_ice.nc
  RSTN=`get_value_in_namelist 'ln_rstart' '&namrun' ${PREFIX}namelist_cfg`
  IS_ICB=`get_value_in_namelist 'ln_icebergs' '&namberg' ${PREFIX}namelist_cfg`
  NIT_RST=${NITENDM1}
  if [ $RSTN == ".true." ]; then
    if [ $NIT_RST -eq 0 ]; then
      ln -s -v ${INPUTDIR}/${PREFIX}${CONFPAR}_restart_00000000.nc ${PREFIX}restart.nc
      ln -s -v ${INPUTDIR}/${PREFIX}${CONFPAR}_restart_ice_00000000.nc ${PREFIX}restart_ice.nc
      if [ $IS_ICB == ".true." ]; then
        ln -s -v ${INPUTDIR}/${PREFIX}${CONFPAR}_restart_icb_00000000.nc ${PREFIX}restart_icb.nc
      fi
    else
      if [ ! -f ${PREFIX}restart_${NIT_RST}.nc ]; then
        echo "Copy ocean restart file from ${STOCKDIR}/restart/nemo_${CONFIG}-${CASE}"
        cp -p ${STOCKDIR}/restart/nemo_${CONFIG}_${CASE}/${PREFIX}restart_${NIT_RST}.nc .
      fi
      ln -s -v ${PREFIX}restart_${NIT_RST}.nc   ${PREFIX}restart.nc
      if [ ! -f ${PREFIX}restart_ice_${NIT_RST}.nc ]; then
        echo "Copy sea ice restart file from ${STOCKDIR}/restart/nemo_${CONFIG}-${CASE}"
        cp -p ${STOCKDIR}/restart/nemo_${CONFIG}_${CASE}/${PREFIX}restart_ice_${NIT_RST}.nc .
      fi
      ln -s -v ${PREFIX}restart_ice_${NIT_RST}.nc   ${PREFIX}restart_ice.nc
      if [ $IS_ICB == ".true." ]; then
        if [ ! -f ${PREFIX}restart_icb_${NIT_RST}.nc ]; then
          echo "Copy iceberg restart file from ${STOCKDIR}/restart/nemo_${CONFIG}-${CASE}"
          if [ -f ${STOCKDIR}/restart/nemo_${CONFIG}_${CASE}/${PREFIX}restart_icb_${NIT_RST}.nc ]; then
            cp -p ${STOCKDIR}/restart/nemo_${CONFIG}_${CASE}/${PREFIX}restart_icb_${NIT_RST}.nc .
            ln -s -v ${PREFIX}restart_icb_${NIT_RST}.nc   ${PREFIX}restart_icb.nc
          elif [ -f ${STOCKDIR}/restart/nemo_${CONFIG}_${CASE}/${PREFIX}restart_icb_${NIT_RST}.tar ]; then
            cp -p ${STOCKDIR}/restart/nemo_${CONFIG}_${CASE}/${PREFIX}restart_icb_${NIT_RST}.tar .
            tar xvf ${PREFIX}restart_icb_${NIT_RST}.tar
            for file in `ls -1 ${PREFIX}${CONFIG}-${CASE}_${NIT_RST}_restart_icb_*.nc` ; do
              f2=`ls -1 $file | sed -e "s/${CONFIG}-${CASE}_${NIT_RST}_//g"`
              mv $file $f2
            done
          fi
        fi
      fi
    fi
  else
    echo "Not in restart mode -> import initial T,S state"
    if [ -z $GLOBAL_SIM ] || [ ! -f ${INPUTDIR}/${PREFIX}istate_TS_${CONFPAR}_${GLOBAL_SIM}_y${Y0}m${M0}.nc ]; then
      ln -s -v ${INPUTDIR}/${PREFIX}istate_TS_${CONFPAR}_y${Y0}m${M0}.nc ${PREFIX}istate_TS.nc
      ln -s -v ${INPUTDIR}/${PREFIX}istate_sea_ice_${CONFPAR}_y${Y0}m${M0}.nc ${PREFIX}istate_sea_ice.nc
    else
      ln -s -v ${INPUTDIR}/${PREFIX}istate_TS_${CONFPAR}_${GLOBAL_SIM}_y${Y0}m${M0}.nc ${PREFIX}istate_TS.nc
      ln -s -v ${INPUTDIR}/${PREFIX}istate_sea_ice_${CONFPAR}_${GLOBAL_SIM}_y${Y0}m${M0}.nc  ${PREFIX}istate_sea_ice.nc
    fi
  fi
  
  #################
  
  echo " "
  echo "Import (links+copy) of input files completed for domain ${iZOOM}."
  echo " "
  
done # iZOOM

echo " "
echo "Launching the long nemo simulation"
echo " "

#=======================================================================================
#=======================================================================================
# 3- Run script
#=======================================================================================
#=======================================================================================

date
echo " "

rm -f app.conf
for iter in `seq 1 $(( SLURM_NTASKS / NB_TASK_PER_NODE ))`; do
   echo "${NB_TASK_NEMO_PER_NODE}-1 ./nemo.exe"        >> app.conf
   echo "${NB_TASK_XIOS_PER_NODE}-1 ./xios_server.exe" >> app.conf
done

if [ $MACHINE == 'occigen' ]; then
  srun --mpi=pmi2 --multi-prog  ./app.conf
elif [ $MACHINE == 'irene' ]; then
  ccc_mprun -f app.conf
else
  echo "ERROR : $MACHINE is an unknown host for this script: adapt srun/ccc_mprun command line."
  echo "   >>>>>>>>>>> STOP !!!!"
  exit
fi

echo " "
date
echo " "

##-- Rebuild mesh_mask, output.abbort and output.init :
rm -f rebuild_nemo.exe rebuild_nemo
ln -s -v ${NEMOdir}/tools/REBUILD_NEMO/BLD/bin/rebuild_nemo.exe
ln -s -v ${NEMOdir}/tools/REBUILD_NEMO/rebuild_nemo
for STUF in mesh_mask bdy_mesh output.abort output.abort_ice output.init output.init_ice; do
  for iZOOM in $(seq 0 ${NZOOM}); do
    if [ $iZOOM -eq 0 ]; then
      PREFIX=''
    else
      PREFIX="${iZOOM}_"
    fi
    if [ -f ${PREFIX}${STUF}_0000.nc ]; then
      NF=`ls -1 ${PREFIX}${STUF}_[0-9][0-9][0-9][0-9].nc |wc -l`
      rebuild_nemo -m -d 1 -x 200 -y 200 -z 1 -t 1 ${PREFIX}$STUF $NF
      rm -f ${PREFIX}${STUF}_[0-9][0-9][0-9][0-9].nc
    fi
  done
done
rm -f nam_rebuild_[0-9][0-9][0-9][0-9][0-9]

##-- export output files:

if [ ! -d OUTPUT_${NRUN} ]; then
  mkdir OUTPUT_${NRUN}
fi

if [ $IS_ISCPL == 'true' ]; then
  if [ ! -f mesh_mask.nc ]; then
    echo '~!@#$%^&* ERROR: mesh_mask.nc not created >>> stop because compulsory in iscpl mode...'
    date
    exit 
  else
    mv mesh_mask.nc OUTPUT_${NRUN}/mesh_mask_${YEAR}${MONTH}${DAY}.nc
  fi
fi

for iZOOM in $(seq 0 ${NZOOM}); do
  if [ $iZOOM -eq 0 ]; then
    PREFIX=''
  else
    PREFIX="${iZOOM}_"
  fi
  mv -f ${PREFIX}${CONFIG}-${CASE}_[1-5][dmy]_*nc OUTPUT_${NRUN}/.
  mv -f ${PREFIX}namelist_ref                     OUTPUT_${NRUN}/${PREFIX}namelist.${NRUN}
  mv -f ${PREFIX}namelist_ice_ref                 OUTPUT_${NRUN}/${PREFIX}namelist_ice.${NRUN}
  mv -f ${PREFIX}ocean.output                     OUTPUT_${NRUN}/${PREFIX}ocean.output.${NRUN}
  rm -f ${PREFIX}namelist_ice_cfg ${PREFIX}namelist_cfg
done

##########################################################
##-- prepare next run if every little thing went all right

NTEST_O=`ls -1 OUTPUT_${NRUN}/${CONFIG}-${CASE}_[1-5][dmy]_*nc |wc -l`
NTEST_R=`ls -1 ${CONFIG}-${CASE}_*_restart_*.nc |wc -l`

FILE_TEST=`ls -1 OUTPUT_${NRUN}/${CONFIG}-${CASE}_*_SBC.nc | tail -1`
NBNAN=`ncdump -v fwfisf ${FILE_TEST} |grep NaN |wc -l`
if [ $NBNAN -gt 0 ]; then
  echo '**************************************************'
  echo '      NaNs were found in the output files         '
  echo '             >>>>>>>>>   stopping here            '
  echo '**************************************************'
  date
  exit
fi 

if [ ${NTEST_O} -gt 0 ] && [ ${NTEST_R} -gt 0 ] && [ $NBNAN -eq 0 ]; then

  if [ $MACHINE == 'occigen' ]; then
    sbatch ./export_nemo_${NRUN}.sh
  elif [ $MACHINE == 'irene' ]; then
    ccc_msub export_nemo_${NRUN}.sh
  fi

  ##-- write last restart time step of mother grid in prod_nemo.db:
  LAST_RESTART_NIT=`ls -1 ${CONFIG}-${CASE}_*_restart_*.nc | tail -1 | sed -e "s/${CONFIG}-${CASE}//g" | cut -d '_' -f2`
  echo " "
  echo "Last restart created at ocean time step ${LAST_RESTART_NIT}"
  echo "  ---> writting this date in prod_nemo.db"
  echo " "
  echo "$LAST_RESTART_NIT" > restart_nit.txt
  ##-- add last restart time step on chidren grids (at the end of last line in prod_nemo.db):
  for iZOOM in $(seq 1 ${NZOOM}); do
    LAST_RESTART_NIT_ZOOM=`ls -1 ${iZOOM}_${CONFIG}-${CASE}_*_restart_*.nc | tail -1 | sed -e "s/${iZOOM}_${CONFIG}-${CASE}//g" | cut -d '_' -f2`
    sed -e "`wc -l prod_nemo.db | cut -d ' ' -f1`s/$/ ${LAST_RESTART_NIT_ZOOM}/g" prod_nemo.db > tmp
    mv tmp prod_nemo.db
    echo " "
    echo "Last restart created for zoom nb $iZOOM at ocean time step ${LAST_RESTART_NIT_ZOOM}"
    echo "  ---> writting this date in prod_nemo.db"
    echo " "
    echo "$LAST_RESTART_NIT_ZOOM" > ${iZOOM}_restart_nit.txt
  done

  echo " "
  date
  echo " "

  ####################################
  ##-- Rebuilding restart files :
  for iZOOM in $(seq 0 ${NZOOM}); do

    if [ $iZOOM -eq 0 ]; then
      PREFIX=''
      LAST_RESTART_NIT_DOM=${LAST_RESTART_NIT}
    else
      PREFIX="${iZOOM}_"
      LAST_RESTART_NIT_DOM=`cat ${iZOOM}_restart_nit.txt`
    fi

    ## ocean and sea-ice restart
    FILEBASE_OCE=`ls -1 ${PREFIX}${CONFIG}-${CASE}_[0-9]???????_restart_0000.nc | sed -e "s/_0000.nc//g"`
    FILEBASE_ICE=`ls -1 ${PREFIX}${CONFIG}-${CASE}_[0-9]???????_restart_ice_0000.nc | sed -e "s/_0000.nc//g"`
    for STUF in $FILEBASE_OCE $FILEBASE_ICE; do
      NF=`ls -1 ${STUF}_[0-9][0-9][0-9][0-9].nc |wc -l`
      rebuild_nemo -m -d 1 -x 200 -y 200 -z 1 -t 1 $STUF $NF
      rm -f ${STUF}_[0-9][0-9][0-9][0-9].nc
    done
    mv ${FILEBASE_OCE}.nc ${PREFIX}restart_${LAST_RESTART_NIT_DOM}.nc
    mv ${FILEBASE_ICE}.nc ${PREFIX}restart_ice_${LAST_RESTART_NIT_DOM}.nc
  
    ## Create an archive for the icb files
    IS_ICB=`get_value_in_namelist 'ln_icebergs' '&namberg' ${PREFIX}namelist_cfg`
    if [ $IS_ICB == ".true." ]; then
      FILEBASE_ICB="${PREFIX}${CONFIG}-${CASE}_${LAST_RESTART_NIT_DOM}_restart_icb.tar"
      tar -cf $FILEBASE_ICB  ${PREFIX}${CONFIG}-${CASE}_${LAST_RESTART_NIT_DOM}_restart_icb*.nc
      if [ ! -f $FILEBASE_ICB ]; then
         echo "Error when creating $FILEBASE_ICB" && exit
      else
         echo "$FILEBASE_ICB ok, deleting individual nectdf files"
         rm -f ${PREFIX}${CONFIG}-${CASE}_${LAST_RESTART_NIT_DOM}_restart_icb*.nc
         mv $FILEBASE_ICB ${PREFIX}restart_icb_${LAST_RESTART_NIT_DOM}.tar
      fi
    fi

  done # iZOOM

  echo " "
  date
  echo " "

  # prepare initial state for following iteration:
  NRUNm1=$NRUN 
  NRUNm2=`expr $NRUN - 1`
  NRUN=`expr $NRUN + 1`
  TMPTMP="${LAST_RESTART_NIT}"
  for iZOOM in $(seq 1 ${NZOOM}); do
    LAST_RESTART_NIT_ZOOM=`cat ${iZOOM}_restart_nit.txt`
    TMPTMP="${TMPTMP} ${LAST_RESTART_NIT_ZOOM}"
  done
  echo "${NRUN} ${YEARf} ${MONTHf} ${DAYf} ${TMPTMP}" >> prod_nemo.db    ## new line

  # clean WORKDIR:
  YEARm2=`expr $YEAR - 2`
  if [ $YEARm2 -lt 1000 ]; then
    YEARm2="0$YEARm2"
  fi
  if [ $IS_BLK_CORE == ".true." ]; then
    for NAMAT in sn_wndi sn_wndj sn_qsr sn_qlw sn_tair sn_humi sn_prec sn_snow sn_slp; do
      ATM_FILE=`grep $NAMAT namelist_cfg | cut -d "'" -f2 | head -1`
      IS_CLIM=`grep $NAMAT namelist_cfg | cut -d ',' -f5 | sed -e "s/ //g" | head -1`
      if [ $IS_CLIM == ".false." ]; then
        rm -f ${ATM_FILE}_y${YEARm2}.nc
        for iZOOM in $(seq 1 $NZOOM); do
           rm -f ${iZOOM}_${ATM_FILE}_y${YEARm2}.nc
        done
      fi
    done
  fi

  rm -f bdyT_tra_y${YEARm2}.nc
  rm -f bdyU_u?d_y${YEARm2}.nc
  rm -f bdyV_u?d_y${YEARm2}.nc
  rm -f bdyT_ice_y${YEARm2}.nc
  rm -f bdyT_ssh_y${YEARm2}.nc


elif [ $NBNAN -gt 0 ]; then

  echo '**************************************************'
  echo '*     NaNs were found in the output files        *'
  echo '*            >>>>>>>>>   stopping here           *'
  echo '**************************************************'
  date
  exit

else

  echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  echo '!@#$%^&* BIG PROBLEM : no output or no restart files created for NEMO !!  >>>>>>> STOP '
  echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  date
  exit

fi

##########################################################
##-- launch next year of the experiment

if [ $MACHINE == 'occigen' ]; then
  sbatch ./run_nemo.sh 
elif [ $MACHINE == 'irene' ]; then
  ccc_msub run_nemo.sh 
fi

echo " "
date
