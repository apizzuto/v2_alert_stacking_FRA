# Stacking v2 alerts using the Fast Response Analysis

This repository is for the application of the fast response analysis to our v2 alert events. The analysis is stacking the results from all of the alert events from the archival v2 alert catalog, beginning in 2011.

For each alert, we perform a point source search which incorporates the millipede scan of our alert events into the likelihood as a penalty term, and we search for additional events in time windows of 1000 s, 2 days, and in all available years of GFUOnline data.

Relies on the [Fast Response Analysis](https://github.com/IceCubeOpenSource/FastResponseAnalysis).

The repository is structured as follows:
* `francis/time_integrated_scripts/`
    - Scripts to search for neutrino sources in the direction of alert events with no temporal signal PDF
* `francis/transient_scripts/`
    - Scripts to search for coincidences within +/- 500 s of each alert event as well as +/- 1 day of alert events
* `francis/universe/`
    - Takes the result of the point source analyses to set constraints on the intrinsic physics properties of populations of neutrino sources
* `francis/icecube_misc/`
    - Information on alert events

The code uses the acronym "FRANCIS", for "**F**ast **R**esponse **A**nalysis for **N**eutrino **C**oincidences from **I**ceCube **S**ignals"