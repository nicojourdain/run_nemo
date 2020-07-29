## Environment to run NEMO simulations and postprocess outputs

Clone the git repository:
```bash
cd ~
git clone XXXXXXXX 
```

Set up correct environment variables:
```bash
cd run_nemo
echo " " >> ~/.bashrc
echo "# Location of templates used for run_nemo scripts" >> ~/.bashrc
echo "export TEMP_NEMO_DIR=`pwd`" >> ~/.bashrc
echo 'export PATH=$TEMP_NEMO_DIR/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
```

# Prepare a simulation

Here we refer to the domain (region, resolution) as "configuration" (e.g., AMUXL12, ORCA025.L75), and to an individual simulation as "case".

```bash
export WORKDIR=$SCRATCHDIR  # choose where you run nemo
mkdir -pv $WORKDIR/run $WORKDIR/output $WORKDIR/restart
```

To prepare a simulation environment from provided configuration templates or from scratch:
```bash
cd $WORKDIR
export NEW_CONFIG='eORCA025.L75'  ## choose your config name
export NEW_CASE='GNJ001'          ## choose your case name
create_sim_env.sh ${NEW_CONFIG} ${NEW_CASE}
```

*NB:* if ``` MY_CONFIG ``` has a corresponding template in ``` $TEMP_NEMO_DIR/template_run ```, it will use the corresponding running scripts and nemalists.

To copy the simulation environment from one of your existiong simulations, do this (if configurations are the same):
```bash
cd $WORKDIR
export SAME_CONFIG='eORCA025.L75'  ## choose your config name
export OLD_CASE='GNJ001''         ## one of your previous simulations
export NEW_CASE='GNJ002'          ## choose your case name
copy_sim_env_from_existing.sh $SAME_CONFIG $OLD_CASE $NEW_CASE
```

If the configurations are different:
```bash
cd $WORKDIR
export OLD_CONFIG='eORCA025.L75'  ## existing config name
export OLD_CASE='GNJ001''         ## one of your previous simulations
export NEW_CONFIG='eORCA025.L75'  ## choose your new config name
export NEW_CASE='GNJ002'          ## choose your new case name
copy_sim_env_from_existing.sh $OLD_CONFIG $OLD_CASE $NEW_CONFIG $NEW_CASE
```

