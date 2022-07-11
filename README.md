# ALMA3
_the plAnetary Love nuMbers cAlculator_

ALMA3 computes loading and tidal Love numbers for a spherically symmetric, radially stratified planet. Both real (time-domain) and complex (frequency-domain) Love numbers can be computed. The planetary structure can include an arbitrary number of layers, and each layer can have a different rheological law. Currently, ALMA3 can model the following linear rheologies: Elastic, Maxwell visco-elastic, Newtonian viscous fluid, Kelvin-Voigt solid, Burgers and Andrade transient rheologies.

## Licensing and attribution

ALMA3 is the result of research work by Daniele Melini, Christelle Saliby and Giorgio Spada. It is distributed free of charge under the BSD 3-clause license, with the hope that it may be useful for research and educational purposes. A copy of the BSD 3-clause license is included in this repository. Please acknowledge the use of ALMA3 by citing the following publication:

Melini, D., Saliby, C. & Spada, G., On computing viscoelastic Love numbers for general planetary models: the ALMA3 code, Geophysical Journal International, 2022, ggac263, https://doi.org/10.1093/gji/ggac263


## Building ALMA3

To build ALMA3, type `cd src` and then `make`. If the build process is successful, a single executable named `alma.exe` will be created in the main directory. ALMA3 requires the [FMLIB multi-precision library](https://dmsmith.lmu.build/) by D.M. Smith, which is included with the ALMA3 package.

## Runnning ALMA3

To run ALMA3, type

`alma.exe <config_file>`

where `<config_file>` is the name of the configuration file containing all the options for the ALMA3 run. Some example configuration files are stored into the `CONFIGS/` folder. Among other parameters, the configuration file contains the name of the data file describing the planetary rheological structure. Some example files are stored into the `MODELS/` folder.

## Contacts

Please refer for any questions about ALMA3 to:

Daniele Melini, Istituto Nazionale di Geofisica e Vulcanologia, Roma, Italy, daniele.melini@ingv.it

Giorgio Spada, University of Bologna, Italy, giorgio.spada@gmail.com



