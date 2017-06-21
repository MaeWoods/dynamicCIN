# dynamicCIN

## Overview
The project is split into simulation and analysis

## Simulation
### Building the simulation code
The code requires GSL to be installed

    cd dynamicCIN/simulation
    make

### Running the simulation code
There are two modes to running simulations: sampling and resimulation. In sampling mode parameters are sampled from the **prior**, in resimulation mode parameters are sampled from a provided file.

#### Running in sampling mode

     ./CINmodel Particles=3 Seed=1 Model=2

The output simulations are written to a csv file *results-SX-PY-MZ.dat*.
The sampled parameters are written to a tab delimited file *parameters-SX-PY-MZ.dat*.

#### Running in resimulation mode
This is automatically selected when a parameter file is provided

    ./CINmodel Particles=3 Seed=1 Model=2 Resim=posterior_pars.txt 

The output simulations are written to a csv file *results-resim-SX-PY-MZ.dat*. Parameters are not written out in the mode.

#### Example
To sample ten particles from the prior, then resample and resimulate. In a shell do the following

    ./CINmodel Particles=10 Seed=1 Model=2
    cat parameters-S1-P* > posterior_pars.txt
    ./CINmodel Particles=10 Seed=1 Model=2 Resim=posterior_pars.txt 