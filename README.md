# Stacking v2 alerts using the Fast Response Analysis

This repository is for the application of the fast response analysis to our v2 alert events. The analysis is stacking the results from all of the alert events from the archival v2 alert catalog, beginning in 2011.

For each alert, we perform a point source search which incorporates the millipede scan of our alert events into the likelihood as a penalty term, and we search for additional events in time windows of 1000 s, 2 days, and in all available years of GFUOnline data.

Relies on the [Fast Response Analysis](https://github.com/IceCubeOpenSource/FastResponseAnalysis). In order to calculate limits on populations of sources, you will also need a version of FIRESONG. Installation of the code and all dependencies is described in `tutorial.md`

The repository is structured as follows:
* `francis/time_integrated_scripts/`
    - Scripts to search for neutrino sources in the direction of alert events with no temporal signal PDF
* `francis/transient_scripts/`
    - Scripts to search for coincidences within +/- 500 s of each alert event as well as +/- 1 day of alert events
* `francis/universe/`
    - Takes the result of the point source analyses to set constraints on the intrinsic physics properties of populations of neutrino sources
* `francis/icecube_misc/`
    - Information on alert events
* `submit_scripts_and_plotting/cluster_jobs/`
    - How all of the trials were submitted to the cluster
* `submit_scripts_and_plotting/followup_plots/`
    - After trials are run, scripts to make sensitivity, fit bias plots
* `submit_scripts_and_plotting/firesong_based_plots/`
    - After scripts run that use `francis.universe`, create sensitivities to populations of neutrino sources

For instructions on how to run each script, including the command line arguments, please see each respective folder. 

### Installation

I've provided a `setup.py` file so that the user may either install via pip (`python -m pip install --editable v2_alert_followup`), or you can append the location of the code to your pythonpath to make sure that the imports work properly. Either way you do this, you should then be able to import the scripts from this analysis using e.g. `import francis.time_integrated_scripts`). 

Interested in running the code but not sure where to start? See `tutorial.md` for a complete walkthrough and a discussion of dependencies.


The code uses the acronym "FRANCIS", for "**F**ast **R**esponse **A**nalysis for **N**eutrino **C**oincidences from **I**ceCube **S**ignals"
