# cfht_mp_tracking
A package of tools to program the CFHT PH2 when tracking minor planets.

This package contains three tools:

`recon_parser.py`  retrieves the names of objects that might are possible occultation tracking targest. Gets the list from the TNORecon web pages.

`minor_planet_ephemeris.py` creates a CFHT PH2 Ephemeris Target JSON object that can be used to program the CFHT using thier PH2 API.

`ph2.py`  builds the complete PH2 submission file, currently I meail that to CFHT and they load it.

*Usage :* 

- Get the target names (TNOs that might occult a star in the given time range)
`recon_parser.py "2018-09-01 00:00:00" "2019-03-31 00:00:00" `

- Compute a CFHT ephemeris for all the target(s)
`./minor_planet_ephemeris.py "2018-09-01 00:00:00" "2019-10-01 00:00:00" "2013 UO17"`

- Upload the observation to CFHT using those ET files as input. List ET files for all the targets of interest on the comamand line.
`./create_cfht_observations.py 18BC11 18BQ03 2013_UO17.txt --bearer-token <bearer_token_from_website>`

