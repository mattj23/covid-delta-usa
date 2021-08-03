# Covid-19 USA Delta Variant Model

## What is this?

This repo contains the beginning of a simplified, mechanics based epidemiological model of SARS-CoV-2 spread in the USA built in order to estimate what the range of plausible outcomes of the delta variant might be.  

This model is probably best described as a monte-carlo agent based simulation.  Though it does not yet take advantage of the typical features for which agent based models are prized (heterogenous populations, geographic features, etc), it simulates individual interactions in order to completely separate the model into a theoretical micro-mechanics portion and an empirical portion, which together in aggregate form the macro behavior.

The advantage of doing things this way is that while the computation becomes far more expensive than fitting statistical models, the micro behavior is implemented in an explicit and obvious way, and so can be trivially modified to take into account almost any imaginable behavior.  The empirical portion provides a means to adjust the overall aggregate macro behavior within the constraints of the micro behavior, and encapsulates social and policy factors which affect the system.

This is in contrast to analytical/semi-theoretical models like compartmental ODEs (SEIR, for example) and statistical methods like Bayesian networks or simple curve-fitting which rely heavily on the empirical fitting of parameters.  The approach taken here is not inherently "better" than other modeling methods, it simply trades complexity in model creation for complexity in computation.

## How the model works

The model centers around a simulated population.  The population may be scaled by an integer value such that one individual in the simulated population represents *n* individuals in the population being modeled.  

Each individual in the population has a set of features which are tracked for them as the simulation progresses through time:
* If they are currently infected with **no** SARS-CoV-2 variant, the **delta** variant, or the **alpha** variant which also serves as a catch-all for any other non-delta variant
* What day in the simulation they were infected with the variant they currently have
* What day in the simulation their symptoms manifest, if they are infected
* What day in the simulation they will get tested for their infection
* Whether or not they have completed vaccination, and if so what day the regimen was complete

Within the simulation, time is tracked in integer days, and the simulation steps forward one day at a time.  At each timestep/day these individuals will interact with each other in stochastic ways, governed by various probability distributions, and some of their interactions will spread the virus.  At the end of each timestep/day certain population statistics will be recorded for later analysis.

Because of the non-deterministic nature of the model, the history is run many times so that statistical tools can be used to analyze the spectrum of histories produced.

### Model initialization 

Because this is a dynamic model, its ultimate behavior is determined by its initial conditions.  As such, the process of initializing the model population before turning the simulator loose is critical.  In order to start a simulation from a realistic place and thus obtain useful results I've chosen to use historical data (real or fake) to play the model forward from the beginning of the pandemic up until the day where the simulation will take over.  

This is where the first modeling choice has to be made: how to handle historical infections.  

There are two problems with naïvely using the documented case record.  
1. First, documented cases are based on positive tests, and these tests are lagged from the date of actual infection.  Unfortunately the lag is relevant because people are contagious well before symptoms form, let alone by the time they manage to get a test which ultimately has a positive result.  Making it more complicated, there is no good data or easy means of estimating the testing lag, as test availability varied wildly during the pandemic and there was often little incentive for asymptomatic or mild cases to get tested.
2. Second, there is an enormous amount of evidence that the number of documented cases seriously underestimates the number of people who actually had infections.  Recently (2021 spring/summer) a series of studies have been conducted examining seroprevalence in US regions, and all are consistent with much higher rates of infection, ranging from the ~3.5x by the CDC<sup>[1](#cdccovidburden)</sup> to roughly 5x for the first 6 months of the pandemic by the NIH<sup>[2](#undiagnosed2021)</sup>. This underestimation is relevant because it affects the rate of natural immunity in the population.

As a result, it doesn't make sense to try to use the historical documented cases directly, as it would produce dramatically different model dynamics.  Since the date and type of infection is the core information in the simulated population, using estimated infection data directly makes the most sense.  The website [covidestim.org](https://covidestim.org)<sup>[3](covidestim)</sup> which has joint affiliation with Harvard and Yale's public health schools and Stanford Medicine, uses what they describe as a "Bayesean evidence synthesis model" to reconstruct estimates of the infection history which would produce the observed test, case, and mortality data.<sup>[4](#chitwood2020)</sup>

The [covidestim.org](https://covidestim.org) data can be retrieved in a single CSV from their website. The entire estimation history is recalculated every time the model is run with new data, so past estimates may change slightly every time the updated data is received.  But it is from this *estimated* historical infections, which runs roughly 3x to 4x the official documented case numbers to account for *both* undercounting and testing lag, which is used to initialize the model.

The initialization process works as follows: 
1. We begin just before the first predicted infections in the US, with the entire population reset to a default initial state, containing no infections
2. We advance one day at a time up until the day *before* the simulation is set to start
    1. We look at the difference between the cumulative number of estimated infections on that day and the current infections in the population, adjusted by the population scale ratio (are we simulating every person or 1:n people)
    2. If we have variant data, we look up the variant proportions for that day. We divide the number of infections which need to be added to the population into their scaled portions, and then we move through the population "infecting" the appropriate number of individuals
    3. The infection mechanics will be covered in further detail later, but to be brief they consist of drawing the date of symptom onset from the appropriate random distribution and marking the individual with the appropriate variant. See the section on [infection and disease progression](#infection-and-disease-progression) for further details.
    4. Additionally, we look the vaccination history provided, and determine the difference between the expected number of completed vaccinations and the number currently in the population.  We select individuals at random from the population (so long as they are more than 30 days from an infection, if they were infected) and mark them as having received the vaccine on the current date.

At the end of this process the population is fully initialized and ready for the simulation.  All future events will be determined by simulation events and not the historical record.


###  Model forward simulation

Once initialized, the historical record is abandoned and the simulation takes over, advancing through the timesteps one by one until the end date of the simulation is reached.

*Note: this mechanic currently only simulates contact between a population and itself. In the future I will be implementing a ratio to account for travel between adjacent populations, however I'm not currently sure what the most efficient way to implement this is, nor what the best way to estimate the cross-population contact probabilities will be.*

1. First, we start by iterating through the entire set of contagious individuals in the population. For each one we apply a normalized **contact probability** to a binomial distribution to determine how many other members of the population they come into contact with.  This is effectively the same as testing each contagious individual against every other individual, flipping a coin with the normalized contact probability to decide if the two meet.
2. The contact probability determines whether or not two individuals in the population have a potentially transmissible encounter, after which the probability of transmission is governed only by the [infection mechanics](#infection-mechanics).  When a member of the population draws a positive integer ***n*** from the binomial distribution based on the contact probability, that value represents the number of other individuals they will come in contact with that day.  We then randomly select ***n*** other individuals from the population to simulate contact with.<sup>[5](#contact_footnote)</sup>
3. For each contact, the probability of infection is calculated from the underlying [infection mechanics](#infection-mechanics), which take into account factors like the carrier's viral variant and time since infection, and the exposed's time-dependent natural/vaccine immunities against the carrier's variant.  A random value is then selected from a uniform distribution and tested against the unique transmission probability for the contact. If it succeeds the transmission event has occurred and the exposed individual is infected with the carrier's variant.
4. Finally, we look to see if we have a vaccination record provided for the day and, if we do, we apply the vaccinations randomly to unvaccinated members of the population who are at least 30 days since being infected.

### The role of the normalized contact probability

While most modeling approaches combine the probability of social contact and the complex biological probability of infection into a single R<sub>0</sub> or growth factor which gets fit based on empirical data, the point of this model is to try to isolate the biological, micro-scale infection mechanics from the separate, macro-scale effects of social behavior.  

The reason to go through all of that effort is to be able to implement any arbitrarily complex micro-scale behavior in a simple and obvious way.  To simulate *any* behavior at the micro-scale, all that needs to be done is to implement the code that enacts the behavior forward in time and (if needed) to add any data fields to the simulated individuals necessary to track related state.  The model can then be run forward many times and the aggregate effect observed.  The cost is a significant increase in computational work.  

The role of the normalized contact probability, which is a value that determines the rate of potentially transmissible events, is to provide a tunable parameter into the model to account for the more difficult to quantify social and policy factors.  This is the value to be adjusted to match empirical data, but one whose effects will always be constrained by the micro-scale behavior.  

Ideally, then, when reasonable parameters for the normalized contact probability are found by examining the full pre-delta history of the US pandemic, they can be preserved and the micro-scale infection behavior can be changed to incorporate the delta variant without requiring a similar amount of post-delta data to gain confidence in the model output.


### Infection mechanics

The infection mechanics are the heart of this model, and they fall into two portions.  The first has to do with how an infection occurs and what happens to the individual in the population who becomes infected.  The second has to do with calculating the probability of infection, taking into account the complex effects of heterogenous and changing immunities.
 
#### Infection and disease progression

A March 2020 paper in the New England Journal of Medicine by Li et al.<sup>[7](#li2020)</sup> attempted to quantify the transmission dynamics of the early pandemic as observed in the first 425 confirmed cases in Wuhan.  Among other things, Li identified both a distribution of the serial interval (an epidimiological term for the generation-to-generation infection interval) and a distribution representing the incubation time (time from infection to onset of symptoms). The latter took the form of an exponentiated Weibull distribution (Li 2020, Figure 2A) which can be described by the parameters `a=12.6245844, c=0.57604443, scale=0.64524305` with no shift from zero.

Several papers also attempted to produce curves for the time dynamics of infectiousness, notably He et al<sup>[8](#he2020)</sup> published May 2020 in Nature and followed by a correction (Ashcroft et al 2020)<sup>[9](#ashcroft2020)</sup>. Ashcroft revisited He's data and produced a corrected gamma distribution with paramters given in Table 1 for the estimate of "infectivity" with respect to days from symptom offset.

Between Li and Ashcroft's correction of He, both the incubation distribution and the "infectivity" distribution are approximated and can together yield a workable estimation of transmission dynamics for SARS-CoV-2 in the simulated model.

Thus, for the non-delta variants, when an individual in the simulated population is "infected", the date of infection is recorded as the date of the contact event which produced the transmission.  The date of symptom onset is then determined by pulling a random value from Li's exponentiated Weibull distribution and adding it to the current simulation day, yielding a future date at which time "symptoms" will manifest.  This is performed once when the individual is "infected" and the date is recorded into the individual's record.

When this now-infected virus carrying individual is part of a future potential transmission event with another individual, the carrier's "infectivity" is determined by determining the days passed since the individual developed symptoms (will be negative for pre-symptomatic transmission) and plugged into Ashcroft/He's gamma distribution.  The infectivity in this case is a scaled, dimensionless probability of transmitting an infection whose absolute magnitude is unimportant<sup>[10](#infectivity)</sup> but whose relative value against other dates and other viral variants will have a meaningful effect on the model output.


---

<a name="cdccovidburden">1</a>: CDC. Cases, Data, and Surveillance. Centers for Disease Control and Prevention. Published February 11, 2020. Accessed August 2, 2021. https://www.cdc.gov/coronavirus/2019-ncov/cases-updates/burden.html

<a name="undiagnosed2021">2</a>: Kalish H, Klumpp-Thomas C, Hunsberger S, et al. Undiagnosed SARS-CoV-2 seropositivity during the first 6 months of the COVID-19 pandemic in the United States. Sci Transl Med. 2021;13(601):eabh3826. doi:10.1126/scitranslmed.abh3826

<a name="covidestim">3</a>: covidestim: COVID-19 nowcasting. Accessed August 2, 2021. https://covidestim.org/

<a name="chitwood2020">4</a>: Chitwood MH, Russi M, Gunasekera K, et al. Reconstructing the Course of the COVID-19 Epidemic over 2020 for US States and Counties: Results of a Bayesian Evidence Synthesis Model. Epidemiology; 2020. doi:10.1101/2020.06.17.20133983

<a name="contact_footnote">5</a>: This is a performance optimization, as it removes an O(n<sup>2</sup>) procedure and replaces it with an O(n) one that accomplishes the same result

<a name="probability_footnote">6</a>:  This probability can be thought of as the conditional P(T|C), where T is the transmitted event and C is the contact event defined by the normalized contact probability.

<a name="li2020">7</a>: Li Q, Guan X, Wu P, et al. Early Transmission Dynamics in Wuhan, China, of Novel Coronavirus–Infected Pneumonia. N Engl J Med. 2020;382(13):1199-1207. doi:10.1056/NEJMoa2001316

<a name="he2020">8</a>: He X, Lau EHY, Wu P, et al. Temporal dynamics in viral shedding and transmissibility of COVID-19. Nat Med. 2020;26(5):672-675. doi:10.1038/s41591-020-0869-5

<a name="ashcroft2020">9</a>: Ashcroft P, Huisman JS, Lehtinen S, et al. COVID-19 infectivity profile correction. Swiss Med Wkly. Published online August 5, 2020. doi:10.4414/smw.2020.20336

<a name="infectivity">10</a>: On an absolute scale the values of the infectivity probability aren't important because the model dynamics are governed by the ratio between the normalized contact probability and the area under the infectivity curve.  That is, if the area under the infectivity curve gets smaller, the contact probability will become larger to accomodate the observed historical record, and vice versa.  The significance of the infectivity curve is that we can then, for example, make a rough estimate of the delta variant's dynamics by plugging in a different curve, or even just a scale factor.

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
---

## C++ and Python components

This project requires both the manipulation/preparation/processing/display of historical and simulation data, and the efficient computation of potentially hundreds of millions of individual interactions between simulated members of a population.

As such, this project is split into two components.  One is a set of python modules and scripts for fetching public data, preparing simulation input data, and retrieving and plotting simulation results.  The other is a C++ program written to efficiently do the simulation itself.

Running the simulations requires both components.  

The python components are located in the `data/` folder of this repository, and scripts should be run out of that directory.  The python scripts are the main way of interacting with the the simulator.

The C++ component is located in the `simulator/` folder.  It should be compiled to a binary before attempting to run the python scripts.

### Preparing the C++ component

The C++ simulator is written in C++20 and requires `nlohmann/json` as a dependency.  OpenMP will be a future requirement for multithreading.  CMake is the build system.

To run this on Linux, you will need to install `nlohman/json` on your machine such that CMake can find it.  At that point you can build the release version of project using CMake

For example, to compile on an Ubuntu 20.04 machine, such as an AWS EC2 instance, follow these instructions:

```bash
# Update out-of-date packages and install dependencies
sudo apt update && apt upgrade -y 
sudo apt install cmake build-essentials nlohmann-json3-dev python3-pip

# Clone the repository
git clone https://github.com/mattj23/covid-delta-usa.git

# Enter the repository, use cmake to build the C++ binary in an out-of-source folder
cd covid-delta-usa
cmake -S simulator/ -B build
cmake --build build/ --target delta_sim --config Release
```

At that point the binary should exist in `~/covid-delta-usa/build/delta_sim`

### Preparing the python component

#### Install the contents of requirements.txt

It's recommended to create a virtual environment to run the python code, but it will work without.  Using `pip` install the contents of `requirements.txt` in the `data/` folder.  This was developed using python 3.9, tested on 3.8, and it may run on versions as old as 3.6 (I have not tested this).  

On an Ubuntu 20.04 machine (once again, using an AWS EC2 instance as an example)
```bash
# Navigate to the data directory of the repository cloned above
cd covid-delta-usa/data

# If you intend to run a virtual environment, create and activate one here
python3 -m venv my_virtenv
source my_vertenv/bin/activate

# Now install the requirements
pip3 install -r requirements.txt
```



#### Set up the settings.json file

In order to retrieve the historical data, including vaccinations and positive tests, you will need an API key for [covidactnow.org](https://covidactnow.org). The API key is free and only takes a few minutes to get.

This key, along with any custom settings, must be put into `settings.json`. Alternatively you can create a file named `settings.gitignore.json` which will supersede `settings.json` and is set to be ignored by the `.gitignore` file.

`settings.json` *or* `settings.gitignore.json`
```json
{
  "act_covid_key": "<act covid api key>",
  "default_input_file": "/tmp/input_data.json",
  "default_output_file": "/tmp/output_data.json",
  "binary_path": "../simulator/cmake-build-release/delta_sim",
  "cache_folder": "./cache"
}
```

Additionally, make sure to set the `binary_path` to point at the binary compiled from the C++ source.  It's best to use absolute paths where possible.

#### Retrieving public data

Once the `settings.json` file has the covidactnow.org API key, run the `fetch_data.py` script to retrieve the actual historical data.  This will save a set of files in the cache folder, also specified in `settings.json`, but defaulting to a directory named `cache/` in the python folder.

```bash
python3 fetch_data.py
```

### Running the model

At this point your environment is prepared and you can run the simulations.  An example is prepared in `simulate.py` which can be run directly.  The file contains a simple setup with a few plots, however they can be modified in any way to investigate any part of the system.  Counter-factual scenarios can be tested by altering input data going to the simulator (for instance, to run tests of different rates of vaccine adoption).

If you aren't running locally, make sure to save figures instead of showing them, or consider setting up a Jupyter notebook.
