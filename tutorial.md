# Getting started

This page walks through the installation of the `v2_alert_stacking` code as well as all of the dependencies. For convenience, we provide a tar file containing an environment with all of the dependencies installed. This can be found at `/data/ana/analyses/NuSources/2021_v2_alert_stacking_FRA/v2_alert_followup.RHEL_7_x86_64.tar.gz`.

Navigate to where you would like this environment to live, and run 
```
cp /data/ana/analyses/NuSources/2021_v2_alert_stacking_FRA/v2_alert_followup.RHEL_7_x86_64.tar.gz ./
tar -xvf v2_alert_followup.RHEL_7_x86_64.tar.gz
```

You can then load the environment by running
```
./combo-plus.RHEL_7_x86_64/env-shell.sh
```

This environment contains the following dependencies:
* CVMFS: `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/setup.sh`
* The realtime metaproject, which is can be built as a parasitic build off of icerec (see instructions on the FRA wiki on how to install this)
* Skylab, [v2.11](https://github.com/icecube/skylab/releases/tag/v2.11)
* [FIRESONG](https://github.com/icecube/FIRESONG): Note, you will need to have a version of FIRESONG where the modules are importable. Right now, this is available in a branch of FIRESONG that can be found in [v1.6](https://github.com/icecube/FIRESONG/releases/tag/v1.6) using the pip-installable branch [here](https://github.com/icecube/FIRESONG/tree/feature/pip-installable).
* [Fast Response Analysis, v0.0.1](https://github.com/icecube/FastResponseAnalysis/releases/tag/v0.0.1)
* [PyCondor](https://github.com/jrbourbeau/pycondor)
    - This is a helpful little package that James Bourbeau wrote for interfacing with the cluster. It can be installed via pip: `pip install pycondor`
* Some more "typical" packages: numpy, scipy, pandas, astropy, healpy, seaborn, matplotlib


### Installing the analysis code
Now, navigate to where the analysis code lives, and run
```
pip install --editable /path/to/2021_v2_alert_stacking_FRA
pip install -r /path/to/2021_v2_alert_stacking_FRA/requirements.txt
```

this should install a module called `francis` into your environment. Most of the requirements are redundant, but it doesn't hurt to double check. 

## Running trials
You made it! You should (hopefully) be able to run the code now! Let's start with how you can run a few trials. These scripts should run fine out of box, but might crash if you try to save the file to a path where you don't have access. 

Want to run some trials for the time-integrated analysis?
```
python francis/time_integrated_scripts/steady_background_trials.py --i=140 --verbose --ntrials=50
```
this will run 50 background trials for the alert event with index=140. it will print some output (because of the verbose flag) and it will account for systematics in the millipede LLH treatment

Want to run short timescale followup trials?
```
python francis/transient_scripts/alert_event_background.py --index=255 --deltaT=1000. --ntrials=100
```

Want to run FIRESONG trials that sample populations of sources?
```
python francis/universe/transient_full_2d_ts_sampling.py --n=50 --density=1e-10 --delta_t=1000.
```

The respective directories include details on all of the scripts and their arguments. The ones here were just examples for each part of the analysis (transient, time-integrated, and stacking)

## Submitting all of the relevant jobs
If you wanted to run all of the trials for these analysis, you would need to use the cluster. All of these scripts are in `submit_scripts_and_plotting/cluster_jobs`. The jobs for the `steady` and `transient` need to be run before `luminosity_constraints`, as the `luminosity_constraints` jobs sample from the trials from the first two.