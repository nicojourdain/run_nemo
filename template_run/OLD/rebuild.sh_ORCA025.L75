#!/bin/bash

# NEMO base directory (you may need to adapt this if it does not work as it is):
WHOAMI=`whoami`
NEMOdir=`grep "NEMOdir=" run_nemo.sh | grep -v "#NEMOdir=" | cut -d '=' -f2 | cut -d '#' -f1 |sed -e "s/'//g ; s/ //g ; s/\\\`whoami\\\`/${WHOAMI}/g ; s/\"//g"`
echo "NEMOdir is $NEMOdir"

cat > tmptmpvbn.txt << EOF
Usage: `basename $0` file_to_rebuild.nc  [zoom_number]              
       file_to_rebuild.nc is the root name of the file to rebuild   
       zoom_number (default=0) is the AGRIF zoom nb to process      
       ex.: ./rebuild.sh mesh_mask                                  
            ./rebuild.sh mesh_mask 2                                
            ./rebuild.sh output.abort 1                             
EOF

#####################################################
# A bit of checking :

if [ ! -n "$1" ]; then
  cat tmptmpvbn.txt
  rm -f tmptmpvbn.txt
  exit
fi

if [ $# -eq 1 ]; then
  FTYPE="$1"
elif [ $# -eq 2 ]; then
  if [ $2 -lt 0 ] || [ $2 -gt 10 ]; then
    cat tmptmpvbn.txt
    rm -f tmptmpvbn.txt
    exit
  fi
  FTYPE="${2}_${1}"
else
  cat tmptmpvbn.txt
  rm -f tmptmpvbn.txt
  exit
fi

#####################################################
# Now rebuilding the MPI subdomains into one file :

rm -rf TMPXUXU
mkdir TMPXUXU
cd TMPXUXU

rm -f rebuild_nemo.exe
ln -s ${NEMOdir}/TOOLS/REBUILD_NEMO/BLD/bin/rebuild_nemo.exe

mv ../${FTYPE}_[0-9]???.nc .
NDOMAIN=`ls -1 ${FTYPE}_[0-9]???.nc |wc -l`

cat > nam_rebuild << EOF
  &nam_rebuild
  filebase='${FTYPE}'
  ndomain=${NDOMAIN}
  /
EOF

cat nam_rebuild
echo " "
echo "./rebuild_nemo.exe"
./rebuild_nemo.exe

mv ${FTYPE}.nc ..
cd ..

if [ -f ${FTYPE}.nc ]; then
  rm -rf TMPXUXU
  echo " "
  echo "${FTYPE}.nc [oK]"
else
  echo "~!@#$%^&* ${FTYPE}.nc has not been created"
fi

rm -f tmptmpvbn.txt
