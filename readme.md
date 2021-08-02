# Covid-19 USA Delta Variant Model

## What is this?

This repo contains the beginning of a simplified, mechanics based epidemiological model of SARS-CoV-2 spread in the USA built in order to estimate what the range of plausible outcomes of the delta variant might be.  

This model is probably best described as a monte-carlo agent based simulation.  Though it does not currently take advantage of the typical features for which agent based models are prized (hetrogenous populations, geographic features, etc), it simulates individual interactions in order to completely separate the model into a theoretical micro-mechanics portion and an empirical portion, which together in aggregate form the macro behavior.

The advantage of doing things this way is that while the computation becomes far more expensive than fitting statistical models, the micro behavior is implemented in an explicit and obvious way, and so can be trivially modified to take into account almost any imaginable behavior.  The empirical portion provides a means to adjust the overall aggregate macro behavior within the constraints of the micro behavior, and encapsualtes social and policy factors which affect the system.

This is in contrast to analytical/semi-theoretical models like SEIR and statistical methods like Bayesian networks or simple curve-fitting which rely heavily on the empirical fitting of parameters.  The approach taken here is not inherently "better" than other modeling methods, it simply trades complexity in model creation for complexity in computation.

## C++ and Python components

This project requires both the manipulation/preparation/processing/display of historical and simulation data, and the efficient computation of potentially hundreds of millions of individual interactions between simulated members of a population.

As such, this project is split into two components.  One is a set of python modules and scripts for fetching public data, preparing simulation input data, and retrieving and plotting simulation results.  The other is a C++ program written to efficiently do the simulation itself.

Running the simulations requires both components.  

The python components are located in the `data/` folder of this repository, and scripts should be run out of that directory.  The python scripts are the main way of interacting with the the simulator.

The C++ component is located in the `simulator/` folder.  It should be compiled to a binary before attempting to run the python scripts.

### Preparing the C++ component

The C++ simulator is written in C++20 and requires `nlohmann/json` as a dependency.  OpenMP will be a future requirement for multithreading.  CMake is the build system.

To run this on Linux, you will need to install `nlohman/json` on your machine such that CMake can find it.  At that point you can build the release version of project using CMake

*Detailed instructions on complilation for linux will be coming*

### Preparing the python component

#### Install the contents of requirements.txt

It's recommended to create a virtual environment to run the python code, but it will work without.  Using `pip` install the contents of `requirements.txt` in the `data/` folder.  This was developed using python 3.9, although it may run on versions as old as 3.6 (I have not tested this).  

#### Set up the settings.json file

In order to retrieve the historical data, including vaccinations and positive tests, you will need an API key for [covidactnow.org](https://covidactnow.org). The API key is free and only takes a few minutes to get.

This key, along with any custom settings, must be put into `settings.json`. Alternatively you can create a file named `settings.gitignore.json` which will supercede `settings.json` and is set to be ignored by the `.gitignore` file.

#### Retriving public data

Once the `settings.json` file has the covidactnow.org API key, run the `fetch_data.py` script to retrieve the actual historical data.  This will save a set of files in the cache folder, also specified in `settings.json`, but defaulting to a directory named `cache/` in the python folder.

