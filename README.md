# ALMA3
_the plAnetary Love nuMbers cAlculator_

ALMA3 computes loading and tidal Love numbers for a spherically symmetric, radially stratified planet. Both real (time-domain) and complex (frequency-domain) Love numbers can be computed. The planetary structure can include an arbitrary number of layers, and each layer can have a different rheological law. Currently, ALMA3 can model the following linear rheologies: Elastic, Maxwell visco-elastic, Newtonian viscous fluid, Kelvin-Voigt solid, Burgers transient rheology and Andrade transient rheology.

## Building ALMA3

To build ALMA3, type `cd src` and then `make`. If the build process is successful, a single executable named `alma.exe` will be created in the main directory. ALMA3 requires the [FMLIB multi-precision library](https://dmsmith.lmu.build/) by D.M. Smith, which is included with the ALMA3 package.

## Runnning ALMA3

To run ALMA3, type

`alma.exe <config_file>`

where `<config_file>` is the name of the configuration file containing all the options for the ALMA3 run. Some example configuration files are stored into the `CONFIGS/` folder. Among other parameters, the configuration file contains the name of the data file describing the planetary rheological structure. Some example files are stored into the `rheological_models/` folder.



