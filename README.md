# micemag - MICE Magnetic Field Mapping Code

Code for modelling the magnetic fields of the spectrometer solenoids used on the Muon Ionisation Cooling Experiment (MICE).

## Getting Started

Clone the code from github to obtain the latest stable release. Although this package comes with a `setup.py` script, it does not work for now.

For pressing ceremonial reasons, this code is written in python 2.7. 


### Prerequisites

For micemag to work you will need to install:

  * scipy
  * numpy
  * matplotlib
  * iminuit

All of which are readily available with `pip`.

## Quick Start Guide

The `micemag/data/` directory contains pre-fitted fields and other data that are used to produce the field maps.  I **highly** recommend that these files are left untouched.  Although it must be said that in future commits these will be updated with better fitted data.

The `bin/` directory contains scripts with neat command line interfaces for out-of-the-box functions such as outputting a magnetic field for use with MAUS.

The main script is `make_g4bl_grid.py` which outputs a g4blGrid format for use in MAUS. The script takes command line arguments for the desired spectrometer solenoid and coil currents. For example:

```
$./make_g4bl_grid.py SSU -M1 185.0 -CC 205.7
```
The script can always be invoked with the `-h` or `--help` option to view usage and available options.

This will make a field map of the SSU magnet with M1 set at 185 A and E1-CC-E2 set at 205.7 A.  Any omitted coils will default to 0 A. Specifying the magnet is mandatory, although if no currents are set the script will exit without doing anything.

The fields are outputted in the `micemag/data/MAUS/` directory and the default names they are saved as are MAGNET_DATE.table, e.g `SSU_20181115.table`.  The fields can get pretty large in size so I recommend that you monitor the amount of fields within that directory.

Right now, the main script only outputs contributions from the Fourier-Bessel expansion but the code is being updated and will become a fuller package.

## Built With

* [numpy](http://www.numpy.org/) - Mathematical library.
* [scipy](https://www.scipy.org/) - Special functions and interpolation routines.
* [matplotlib](https://matplotlib.org/) - Used for beautiful plots.
* [iminiuit](https://github.com/scikit-hep/iminuit) - Used for minimization.

## Authors

* **Joe Langlands** - *Current maintainer* - [JoeLanglands](https://github.com/JoeLanglands)
* **Victoria Blackmore** - *Initial contributor*

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Victoria Blackmore for writing all the base code that eventually evolved into this.
