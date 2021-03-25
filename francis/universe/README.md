# Scripts for looking at populations of sources

This directory contains the scripts to run trials and to make plots for seeing how the analysis would respond to various populations of neutrino sources. 

Scripts:
`steady_2d_ts_sampling.py`: Given a density and luminosity function of neutrino sources, calculate the expected TS and binomial p-value distribution. Args:
* `--n`: Number of trials (default 1000)
* `--density`: Local source density in Mpc^-3
* `--LF`: Luminosity function of neutrino sources to pass to FIRESONG
* `--evol`: Evolution of sources to assume. Default to Madau and Dickinson SFH model ("MD2014SFR)
* `--manual_lumi`: If 0.0, then find the luminosity that saturates the diffuse flux based on the population. Else, use the given luminosity (in same units as FIRESONG)

`transient_full_2d_ts_sampling.py`: Given a density and luminosity function of neutrino sources, calculate the expected TS and binomial p-value distribution. Same as `steady_2d_ts_sampling`, but for transient sources. Args:
* `--n`: Number of trials (default 1000)
* `--density`: Local source density in Mpc^-3
* `--LF`: Luminosity function of neutrino sources to pass to FIRESONG
* `--evol`: Evolution of sources to assume. Default to Madau and Dickinson SFH model ("MD2014SFR)
* `--manual_lumi`: If 0.0, then find the luminosity that saturates the diffuse flux based on the population. Else, use the given luminosity (in same units as FIRESONG)
* `--delta_t`: Timescale of the analysis to use. Either 1000. or 172800

The other `.py` files in this directory are not meant to be run, but they facilitate the running of the scripts above, as well as handle the plotting of relevant analysis features. I provide a quick summary of what they do, below:

`transient_universe.py`: based off of an assumed density/luminosity of a source population, find alert events. For every signal alert event, find those that will be accompanied by additional events in the GFU sample

`universe_analysis.py`: From the alerts found in `transient_universe`, sample the relevant TS distributions, to find the analysis response to a given population of sources.

`universe_plotter.py`: From the analyses performed using `universe_analysis`, calculate sensitivities to populations of sources.