# Environment to run NEMO simulations and postprocess outputs

_Contributors:_
* Nicolas Jourdain (IGE-CNRS)

_Tests:_
* Has been tested for NEMO-3.6 (tag "r3.6") and NEMO-4.2-RC (tag "r4.2.RC").

_Known caveats:_
* Currently only works on occigen (CINES) and Irene-Rome (TGCC).


To clone the git repository (the standard version is made for NEMO4.2.0):
```bash
cd ~ # or anywhere else
git clone git@github.com:nicojourdain/run_nemo.git # if you use github with SSH key
git clone https://github.com/nicojourdain/run_nemo.git # otherwise
```
If instead you use NEMO-3.6 (tag 'r3.6') and NEMO-4.2-RC (tag 'r4.2.RC'), clone as follows:
```bash
cd ~ # or anywhere else
BRANCH='r3.6' # or 'r4.2.RC'
git clone --branch $BRANCH git@github.com:nicojourdain/run_nemo.git # if you use github with SSH key
git clone --branch $BRANCH https://github.com/nicojourdain/run_nemo.git # otherwise
```

Set up correct environment variables (do it once):
```bash
cd run_nemo
echo " " >> ~/.bashrc
echo "# Location of templates used for run_nemo scripts" >> ~/.bashrc
echo "export TEMP_NEMO_DIR=`pwd`" >> ~/.bashrc
echo 'export PATH=$TEMP_NEMO_DIR/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
```

## Prepare the environment for a simulation

Here we refer to the domain (region, resolution) as "configuration" (e.g., AMUXL025.L75, AMUXL12, ORCA025.L75), and to an individual simulation as "case".

```bash
export SCRATCHDIR=${GEN6035_CCCSCRATCHDIR}  # choose where you run nemo
mkdir -pv $SCRATCHDIR/run $SCRATCHDIR/output $SCRATCHDIR/restart
```

To prepare a simulation environment from provided configuration templates or from scratch:
```bash
cd $SCRATCHDIR
export NEW_CONFIG='AMUXL025.L75'  ## choose your config name
export NEW_CASE='TEST01'          ## choose your case name
create_sim_env.sh ${NEW_CONFIG} ${NEW_CASE}
```

*NB:* if ``` MY_CONFIG ``` has a corresponding template in ``` $TEMP_NEMO_DIR/template_run ```, it will use the corresponding running scripts and namelists.

To copy the simulation environment from one of your existiong simulations, do this (if configurations are the same):
```bash
cd $SCRATCHDIR
export SAME_CONFIG='eORCA025.L75'  ## choose your config name
export OLD_CASE='GNJ001''         ## one of your previous simulations
copy_sim_env_from_existing.sh $SAME_CONFIG $OLD_CASE $NEW_CASE
```

If the configurations are different:
```bash
cd $SCRATCHDIR
export OLD_CONFIG='eORCA025.L121'  ## existing config name
export OLD_CASE='GNJ001''          ## one of your previous simulations
copy_sim_env_from_existing.sh $OLD_CONFIG $OLD_CASE $NEW_CONFIG $NEW_CASE
```

## Run a NEMO simulation

To prepare the input files, see [https://github.com/nicojourdain/BUILD_CONFIG_NEMO](https://github.com/nicojourdain/BUILD_CONFIG_NEMO). It is recommended to store the input files in ``` $SCRATCHDIR/input/nemo_${NEW_CONFIG} ```. 

Then, to run NEMO:

```bash
cd $SCRATCHDIR/run/nemo_${NEW_CONFIG}_${NEW_CASE}
vi run_nemo.sh # edit User's choices
vi namelist_nemo_GENERIC_${NEW_CONFIG} # make your choices
vi namelist_ice_nemo_GENERIC_${NEW_CONFIG} # make your choices
vi iodef_monthly_daily.xml # edit output variables (if you choose BY_MONTH=1)
vi iodef_daily.xml # edit output variables (if you choose BY_MONTH=0)
```

Then, to launch the simulation:
```bash
sbatch run_nemo.sh
```

The simulation is resubmitted every NDAYS until it reaches ``` YEAR_MAX ``` or a total number of ```NRUN_MAX ``` resubmissions. The progress of the overall simulation is kept in ``` prod_nemo.db ```.

After each submission, outputs are temporarilly stored in ``` $SCRATCHDIR/run/nemo_${NEW_CONFIG}_${NEW_CASE}/OUTPUT_xx ``` before being treated by ``` compress_nemo_xx.sh ```. Then, the outputs and the used namelists are stored in ``` $SCRATCHDIR/output/nemo_${NEW_CONFIG}_${NEW_CASE} ``` and the restart files in ``` $SCRATCHDIR/restart/nemo_${NEW_CONFIG}_${NEW_CASE} ```.

## Postprocess the outputs

```bash
cd $SCRATCHDIR/output/nemo_${NEW_CONFIG}_${NEW_CASE}
vi postprocess_nemo_with_monthly_io.sh # edit (if you used BY_MONTH=1)
sbatch postprocess_nemo_with_monthly_io.sh
```

