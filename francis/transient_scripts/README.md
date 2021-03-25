# Transient followup scripts

This directory contains the scripts to run trials for the short timescale followups for each alert event included in the v2 alert catalog. 

Scripts:
`alert_event_background.py`: Runs background trials. Args:
* `--ntrials`: Number of trials (default 1000)
* `--index`: Alert event index, identifier for alert events ordered by date (from 0 to 276)
* `--smear`: Raise this flag if you want to account for systematics when converting the millipede LLH map to a spatial prior
* `--deltaT`: Time window in seconds (use either 1000 or 172800)

`alert_event_fits.py`: Runs round-trip injection tests. Args:
`alert_event_background.py`: Runs background trials. Args:
* `--ntrials`: Number of trials per signal strength (default 100)
* `--index`: Alert event index, identifier for alert events ordered by date (from 0 to 276)
* `--smear`: Raise this flag if you want to account for systematics when converting the millipede LLH map to a spatial prior
* `--deltaT`: Time window in seconds (use either 1000 or 172800)

`alert_event_sensitivity.py`: Runs signal injection trials for sensitivity. Args:
* `--ntrials`: Number of trials per signal strength (default 100)
* `--index`: Alert event index, identifier for alert events ordered by date (from 0 to 276)
* `--smear`: Raise this flag if you want to account for systematics when converting the millipede LLH map to a spatial prior
* `--deltaT`: Time window in seconds (use either 1000 or 172800)

`alert_event_unblind.py`: Looks at actual data for one alert event. Args:
* `--index`: Alert event index, identifier for alert events ordered by date (from 0 to 276)
* `--smear`: Raise this flag if you want to account for systematics when converting the millipede LLH map to a spatial prior
* `--deltaT`: Time window in seconds (use either 1000 or 172800)

`alert_event_run.py`: Similar to the unblinding script, but this is the script we would run in realtime when new alerts come in. Syntax for arguments similar to the FRA