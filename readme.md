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

### Running the model

See `simulate.py` for an example of running the model.  More information to come.

## Sample model output

The model output is the accumulated historical trajectories of each simulation.  A set of helper classes exist to parse the data and make it easy to plot and analyze.  As an example, here's a run for Florida based on the contents of `simulate.py`.

Keep in mind this is not a hard prediction for FL's future.  The model at this point (2021-07-30) makes many assumptions which may or may not reflect reality.

* The model uses [covestim.org](https://covestim.org)'s estimates of true infections based on a Bayesian inference model, which currently estimates that there have been several times the amount documented of covid cases in the USA.  Estimates of actual cases in the US range from >3x (CDC) documented cases to ~5x (NIH) based on random seropositivity testing.  This is an *optimistic* assumption, because it implies there is a higher amount of natural immunity in the US population than the raw case numbers would suggest.
* The model assumes that a case of the alpha variant confers 100% immunity to all SARS-CoV-2 variants forever.  This is an unrealistic and overly optimistic assumption, which we know is not true. It's currently the next thing on the list to be updated.
* The model assumes that all vaccinations were Pfizer/Moderna and that there is no difference in immunity between the alpha and delta variants.  We already know from a study in the UK<sup>[1](#lopezbernal2021)</sup> this assumption is unrealistic, and the Israeli Ministry of Health data raises questions about time effects on vaccine immunity against delta.  This will also be addressed soon.
* The model does not currently try to predict and thus take into account vaccinations which occur after the historical data from covidactnow.org. If one were to have predicted numbers for the future they are trivial to feed into the simulation, but the image below does not contain an attempt to do so.  Know that that is a *pessimistic* assumption.
* Currently there isn't a good single estimate of the difference in infectivity between the delta variant and prior variants.  This data can be estimated from the observed relative spread rates of the delta and alpha variants in the USA between May and July if reasonable estimates for both natural and vaccine immunities to delta can be made.  I intend to do that soon, but it has not yet been performed.  In lieu of this the image below uses a scale factor of 3.8 on the infectivity curve of the alpha variant<sup>[2](#ashcroft2020)</sup> because it allows the empirical contact probability to remain relatively stable while still accounting for late July 2021 infections.
* This doesn't take into account changes in social contact, which may be an *optimistic* assumption depending on schools and colleges going back into session resulting in more contact than the model assumes, or a *pessimistic* one if people begin social distancing, isolation, and other self-protective activities would would result in less contact than the model assumes.

![FL Simulation Test](https://github.com/mattj23/covid-delta-usa/blob/main/images/fl_test_plot.png?raw=true)

*The methodology used to create this graph was to first look back at late May 2021, before widespread incidence of the delta variant, and adjust the contact probability until it roughly matched case growth at that time.  Then, jumping forward to mid July 2021, leaving the contact probability the same as early summer and adjusting the delta variant infectivity ratio until the infection uptick at the end of July was matched. From there allow the simulation to run out until November, when the system has saturated.*

---


<a name="lopezbernal2021">1</a>: 1. Lopez Bernal J, Andrews N, Gower C, et al. Effectiveness of Covid-19 Vaccines against the B.1.617.2 (Delta) Variant. N Engl J Med. Published online July 21, 2021:NEJMoa2108891. doi:10.1056/NEJMoa2108891

<a name="ashcroft2020">2</a>: Ashcroft P, Huisman JS, Lehtinen S, et al. COVID-19 infectivity profile correction. Swiss Med Wkly. Published online August 5, 2020. doi:10.4414/smw.2020.20336
