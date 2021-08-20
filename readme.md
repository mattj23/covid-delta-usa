# Covid-19 USA Delta Variant Model

## What is this?

This repo contains the beginning of a simplified, mechanics based epidemiological model of SARS-CoV-2 spread in the USA built in order to estimate what the range of plausible outcomes of the delta variant might be.  

This model is probably best described as a monte-carlo agent based simulation.  Though it does not yet take advantage of the typical features for which agent based models are prized (heterogenous populations, geographic features, etc), it simulates individual interactions in order to completely separate the model into a micro-scale portion which behaves according to game-like mechanics, and a macro-scale empirical portion, which together in aggregate form the overall system behavior.

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

Because this is a dynamic model, its behavior is determined by its initial conditions.  As such, the process of initializing the model population before turning the simulator loose is critical.  In order to start a simulation from a realistic place and thus obtain useful results I've chosen to use historical data (real or fake) to play the model forward from the beginning of the pandemic up until the day where the simulation will take over.  

This is where the first modeling choice has to be made: how to handle historical infections.  

There are two problems with naïvely using the documented case record.  
1. First, documented cases are based on positive tests, and these tests are lagged from the date of actual infection.  Unfortunately the lag is relevant because people are contagious well before symptoms form, let alone by the time they manage to get a test which ultimately has a positive result.  Making it more complicated, there is no good data or easy means of estimating the testing lag, as test availability varied wildly during the pandemic and there was often little incentive for asymptomatic or mild cases to get tested.
2. Second, there is an enormous amount of evidence that the number of documented cases seriously underestimates the number of people who actually had infections.  Recently (2021 spring/summer) a series of studies have been conducted examining seroprevalence in US regions, and all are consistent with much higher rates of infection, ranging from the ~3.5x by the CDC<sup>[1](#cdccovidburden)</sup> to roughly 5x for the first 6 months of the pandemic by the NIH<sup>[2](#undiagnosed2021)</sup>. This underestimation is relevant because it affects the rate of natural immunity in the population.

As a result, it doesn't make sense to try to use the historical documented cases directly, as it would produce dramatically different model dynamics.  Since the date and type of infection is the core information in the simulated population, using estimated infection data directly makes the most sense.  The website [covidestim.org](https://covidestim.org)<sup>[3](covidestim)</sup> which has joint affiliation with Harvard and Yale's public health schools and Stanford Medicine, uses what they describe as a "Bayesean evidence synthesis model" to reconstruct estimates of the infection history which would produce the observed test, case, and mortality data.<sup>[4](#chitwood2020)</sup>

The [covidestim.org](https://covidestim.org) data can be retrieved in a single CSV from their website. The entire estimation history is recalculated every time the model is run with new data, so past estimates may change slightly every time the updated data is received.  But it is from this *estimated* historical infections, which runs roughly 3x to 4x the official documented case numbers to account for *both* undercounting and testing lag, which is used to initialize the model.

The initialization process works as follows: 
1. We begin just before the first predicted infections in the US, with the entire population reset to a default initial state, containing no infections
2. We advance one day at a time up until the day *before* the simulation is set to start
    1. We look at the difference between the cumulative number of estimated infections on that day and the current infections in the population, adjusted by the population scale ratio (are we simulating every person or 1:n people?)
    2. If we have variant data, we look up the variant proportions for that day. We divide the number of infections which need to be added to the population into their scaled portions, and then we move through the population "infecting" the appropriate number of individuals with each variant
    3. The infection mechanics will be covered in further detail later, but to be brief there is both natural and vaccine-induced immunities that must be accounted for even in the population initalization. The immunity odds ratios are respected during this phase, making even the population initialization somewhat non-determinisitic. 
    4. Once an individual is selected for infection, the date of their symptom onset is drawn from the appropriate random distribution and the individual is marked with the appropriate variant. A natural immunity scalar is selected for future re-exposures.  See the section on [infection and disease progression](#infection-and-disease-progression) for further details.
    4. Additionally, we look the vaccination history provided, and determine the difference between the expected number of completed vaccinations and the number currently in the population.  We select individuals at random from the population (so long as they are more than 30 days from an infection, if they were infected) and mark them as having received the vaccine on the current date.  They are assigned a random vaccine immunity scalar at this point as well.

At the end of this process the population is fully initialized and ready for the simulation.  All future events will be determined by simulation events and not the historical record.

A copy of the initialized population is maintained so that repeated runs of the simulation do not need to repeat the initialization process.


###  Model forward simulation

Once initialized, the historical record is abandoned and the simulation takes over, advancing through the timesteps one by one until the end date of the simulation is reached.

*Note: this mechanic currently only simulates contact between a population and itself. In the future I will be implementing a ratio to account for travel between adjacent populations, however I'm not currently sure what the most efficient way to implement this is, nor what the best way to estimate the cross-population contact probabilities will be.*

1. First, we start by iterating through the entire set of contagious individuals in the population. For each one we apply a normalized **contact probability** to a binomial distribution to determine how many other members of the population they come into contact with.  This is effectively the same as testing each contagious individual against every other individual, flipping a coin with the normalized contact probability to decide if the two meet.
2. The contact probability determines whether or not two individuals in the population have a potentially transmissible encounter, after which the probability of transmission is governed only by the [infection mechanics](#infection-mechanics).  When a member of the population draws a positive integer ***n*** from the binomial distribution based on the contact probability, that value represents the number of other individuals they will come in contact with that day.  We then randomly select ***n*** other individuals from the population to simulate contact with.<sup>[5](#contact_footnote)</sup>
3. For each contact, the probability of infection is calculated from the underlying [infection mechanics](#infection-mechanics), which take into account factors like the carrier's viral variant and time since infection, and the exposed's time-dependent natural/vaccine immunities against the carrier's variant.  A random value is then selected from a uniform distribution and tested against the unique transmission probability for the contact. If it succeeds the transmission event has occurred and the exposed individual is infected with the carrier's variant.
4. Finally, we look to see if we have a vaccination record provided for the day and, if we do, we apply the vaccinations randomly to unvaccinated members of the population who are at least 30 days since being infected.

### The role of the normalized contact probability

While most modeling approaches combine the probability of social contact and the complex biological probability of infection into a single R<sub>0</sub> or growth factor which gets fit based on empirical data, the point of this model is to try to isolate the biological, micro-scale infection mechanics from the separate, macro-scale effects of social behavior.  

The reason to go through all of that effort is to be able to implement any arbitrarily complex micro-scale behavior in a simple and obvious way.  To simulate *any* behavior at the micro-scale, all that needs to be done is to implement the code that enacts the behavior forward in time and (if needed) to add any data fields to the simulated individuals necessary to track related state.  The model can then be run forward many times and the aggregate effect observed.  The cost of this profound simplicity is a significant increase in CPU time (*this tradeoff is somewhat mitigated by writing extremely performance-conscious C++ code, but these simulations are still computationally expensive compared to curve fitting or solving PDEs*).

The role of the normalized contact probability, which is a value that determines the rate of potentially transmissible events, is to provide a tunable parameter into the model to account for the more difficult to quantify social and policy factors.  This is the value to be adjusted to match empirical data, but one whose effects will always be constrained by the micro-scale behavior.  

Ideally, then, when reasonable parameters for the normalized contact probability are found by examining the full pre-delta history of the US pandemic, they can be preserved and the micro-scale infection behavior can be changed to incorporate the delta variant without requiring a similar amount of post-delta data to gain confidence in the model output.


### Infection mechanics

The infection mechanics are the heart of this model, and they fall into two portions.  The first has to do with how an infection occurs and what happens to the individual in the population who becomes infected.  The second has to do with calculating the probability of infection, taking into account the complex effects of heterogenous and changing immunities.

***Note**: the point of this model is to be able to change any micro-scale behavior arbitrarily.  The following section describes the various probability distributions and scalar curves *I* chose to use in the default infection mechanics, but *any* parameter or curve can be trivially swapped for something different.*
 
#### Infection and disease progression

A March 2020 paper in the New England Journal of Medicine by Li et al.<sup>[7](#li2020)</sup> attempted to quantify the transmission dynamics of the early pandemic as observed in the first 425 confirmed cases in Wuhan.  Among other things, Li identified both a distribution of the serial interval (an epidimiological term for the generation-to-generation infection interval) and a distribution representing the incubation time (time from infection to onset of symptoms). The latter took the form of an exponentiated Weibull distribution (Li 2020, Figure 2A) which can be described by the parameters `a=12.6245844, c=0.57604443, scale=0.64524305` with no shift from zero.

Several papers also attempted to produce curves for the time dynamics of infectiousness, notably He et al<sup>[8](#he2020)</sup> published May 2020 in Nature and followed by a correction (Ashcroft et al 2020)<sup>[9](#ashcroft2020)</sup>. Ashcroft revisited He's data and produced a corrected gamma distribution with paramters given in Table 1 for the estimate of "infectivity" with respect to days from symptom offset.

Between Li and Ashcroft's correction of He, both the incubation distribution and the "infectivity" distribution are approximated and can together yield a workable estimation of transmission dynamics for SARS-CoV-2 in the simulated model.

| ![Infection/Infectivity Distributions](https://github.com/mattj23/covid-delta-usa/blob/main/images/infection_distributions.png?raw=true)|
|:--:|
| *On the left is the current exponentiated weibull distributions used for generating the day of symptom onset for the two model variants, and on the right are the relative infectivity curves used to determine how infectious a simulated individual is on any arbitrary day in the simulation.  These are the default model parameters, but are easily replaced.*|

Thus, for the non-delta variants, when an individual in the simulated population is "infected", the date of infection is recorded as the date of the contact event which produced the transmission.  The date of symptom onset is then determined by pulling a random value from Li's exponentiated Weibull distribution and adding it to the current simulation day, yielding a future date at which time "symptoms" will manifest.  This is performed once when the individual is "infected" and the date is recorded into the individual's record.

When this now-infected virus carrying individual is part of a future potential transmission event with another individual, the carrier's "infectivity" is determined by determining the days passed since the individual developed symptoms (will be negative for pre-symptomatic transmission) and plugged into Ashcroft/He's gamma distribution.  The infectivity in this case is a scaled, dimensionless probability of transmitting an infection whose absolute magnitude is unimportant<sup>[10](#infectivity)</sup> but whose relative value against other dates and other viral variants will have a meaningful effect on the model output.

For the delta variant, a very recent July 23, 2021 preprint<sup>[11](#li2021)</sup> studying an outbreak in Guangzhou noted that relative viral loads were 1260x higher on the day of first detection than they were for older variants.  Additionally, the day of the first positive PCR+ test was 4 days for the delta variant vs 6 days for the 2020 outbreak.  Thus, for a quick and dirty estimation, this model estimates the incubation time of the delta variant by drawing from Li's exponentiated Weibull and then scaling the result by 2/3rds.

The >1000-fold increase in viral load does not correspond with a 1000x increase in infectivity, but rough estimates of R<sub>0</sub> values for traditional models are placing delta at 8 to 9<sup>[12](#cdcleak)</sup>, while the pre-alpha variants were around 3 and the alpha variant is potentially 4 to 5.  As such it's possible to make a rough estimate of delta's infectivity by scaling the output of Ashcroft/He's gamma distribution between 2 and 3. More detailed estimation is covered in the section below devoted to discovering the [relative infectivity of the delta variant.](#relative-infectivity-of-delta-variant)

Additionally, other probability distributions can and will be added here to track disease progression; this will allow the recording of aggregate population statistics for things like disease severity, hospitalizations, and mortality.

#### Immunity and probability of infection

On the disease carrier side, the "infectivity" is computed as described above.  If the exposed individual has an immunity, however, the probability check against the carrier's infectivity is pre-empted and the transmission event is not completed.

*Note: by pre-empting the infectivity check before it happens based on a binary value of immunity, the micro-mechanics explicitly do not contain any mechanism by which a person's immunity reacts differently to various levels of infective potency of the carrier. This could be implemented in the future if such a mechanism could be characterized.*

Immunity is a complex subject even in this simplified model.  Unfortunately, there are no good direct empirical measurements of immunity isolated from social and other factors.  The only measures of immunity we have are efficacy values which represent a reduced incidence ratio between groups.  Because the system being studied is **not** under steady state conditions, efficacy values are products of the exact time period in which they were measured and do not necessarily represent direct underlying probabilities<sup>[13](#efficacy_footnote)</sup>.

To prevent an immunity mechanic that wanes by nature of being repeatedly testing, it's important to not use a "saving throw" mechanism for pre-empting the carrier's infectivity probability check.  Rather, the stochastic portion of an individual's immunity should be resolved once, so that the simulation isn't doing the equivalent of asking over and over again until it gets an answer it likes.  

However, it's also important to be able to account for immunity that changes with time, since it's known<sup>[14](#lopezbernal2021)</sup> that mRNA vaccine immunity changes between the first and second dose and there's evidence<sup>[12](#cdcleak)</sup> that natural immunity against the delta variant may wane beyond 180 days.

To account for both of these seemingly contradictory requirements, the following method is used:
1. A baseline curve is produced for every time-dependant immunity response. This curve is the complement of the odds-ratio based efficacy at any point in time. That is, if the observed population efficacy is 30% at time ***t<sub>0</sub>***, the curve at time ***t<sub>0</sub>*** has a value of 0.7.
2. When an individual joins an immunity category, that individual is assigned a single random value pulled from a uniform distribution from 0.0 to 1.0, called their *immunity scalar*.  That value does not change again.
3. At any given time ***t***, the individual is immune if their *immunity scalar* is greater than the baseline curve at time ***t***.

Assuming that the uniform distribution generates random samples correctly, at any given time 30% of the population who are at ***t<sub>0</sub>*** in their immunity progression will be immune and the other 70% will not.  If the baseline immunity curve rises to 88% at ***t<sub>1</sub>***, any individuals who were already immune will stay immune, some individuals who were not immune will become immune, and some individuals who were not immune will remain without immunity...but overall the proportion of immune individuals will have risen from 30% to 88% without any recalcuation of the *immunity scalar*.

| ![Immunity curves](https://github.com/mattj23/covid-delta-usa/blob/main/images/immunity_curves.png?raw=true)|
|:--:|
| *The natural immunity curves in the left figure are current model defaults. The vaccine efficacy curve from Thomas et al. 2021 is the default used to simulate immunity against the alpha variant. Either the Lopez Bernal et al. curve or the version adjusted with the Israeli Ministry of Health immunity decay can be used in the model to simulate vaccine immunity vs delta, and as expected they produce very different results.*|

For estimating the natural immunity conferred by alpha against alpha, the clearest data I found was Hall et al.<sup>[15](#hall2021)</sup> which tracked a cohort of healthcare providers in the UK which had seropositivity tests done from early in the pandemic. From that we can make at least a rough, usable estimate of the efficacy of natural immunity.

Estimates of infection-acquired immunity against the delta variant are harder to come by.  Public Health England's July 9th, 2021 Technical Briefing 19 offered an estimate of an adjusted odds ratio of reinfection by delta against reinfection by alpha.  From this, modifications can be made to the curve from Hall to produce a similar curve.  Unfortunately, there is some contradiction between the two sources in the time estimates given as to when the odds of reinfection become significant, so this portion of the model could use some attention and careful re-reading of both documents.

Estimates of the Pfizer vaccine efficacy (used as a stand-in for both MRNA based vaccines, which are the dominant vaccines administered in the USA) against the alpha variant and ancestral strains come from Thomas et al. 2021<sup>[16](#thomas2021)</sup>, which is currently a preprint but contains a figure showing efficacy development over time.

The efficacy of the vaccines against the delta variant, however, is much less clear.

There are two sources of data currently (July/August 2021) circulating.  The first is Lopez Bernal et al.<sup>[14](#lopezbernal2021)</sup>, a test-negative case-control study done in the UK.  They calculated overall efficacy of the BNT162b2 vaccine (Pfizer/BioNTech) as 93.7% against the alpha variant and 88.0% against the delta variant.

The other source is a recently released document by the Israeli Ministry of Health<sup>[16](#israelimoh2021)</sup> which implies significant waning/decaying immunity against the delta variant. This document has received a fair amount of criticism in light of the optimistic UK data, as it implies efficacies waning to as low as 16% in the earliest vaccinated cohort.  The Israeli data is based on much smaller case counts than the UK data.

| ![Israeli MOH](https://github.com/mattj23/covid-delta-usa/blob/main/images/israeli_moh_2021-06-20.png?raw=true)|
|:--:|
| *This is from the final slide in the Israeli MOH document, showing decreasing efficacy in earlier cohorts.*|

However, I think the Israeli data is not worth dismissing for the following reasons:
- The UK study is based on vaccinations and positive PCR tests that were performed up to May 16, 2021, and only included symptomatic positive cases.  
- The UK's [2nd dose vacciations](https://coronavirus.data.gov.uk/details/vaccinations) were less than 1M before the start of March, reached 5M by the end of March, 15M by the end of April, and 20M by May 16. That is to say that 50% of the UK's Pfizer-vaccinated population as of the end of the study received their second dose in April and another 25% in May.
- To count as having had 2 doses, a positive case had to have symptoms occuring 14 days or more after the receipt of the second dose, which effectively means that all of the study's 2-dose Pfizer cases had been vaccinated in April or earlier and had to have reported symptoms and tested positive by May 16.
- Thus for the UK study, of the 15M who were fully vaccinated in time to be counted as 2-dose Pfizer cases by the end of the study, 67% had been fully vaccinated for less than 1.5 months, 26% more for less than 2.5 months, and only 6% of them had been vaccinated more than 2.5 months. They were all counted in the same bin when efficacy was calculated.
- Unlike the UK study which stopped collecting data in May, the Israeli data contains cases which occurred through the middle of July
- Israel had fully vaccinated 1.8M by the end of January, 3.3M by the end of February, 4.8M by the end of March, and then reached 5M around the end of April. By the time the study ended in mid July, they had reached 5.2M full vaccinations. That is to say that of the people vaccinated by the end of their study, 34% were vaccinated for more than 5.5 months, another 29% for 4.5 months, another 29% 3.5 months, and only 8% for less than 3.5 months.
- Israel's cohort fully vaccinated in April who had *symptomatic* infections yielded a 79% vaccine efficacy. They were vaccinated at the same time as 67% of the UK cohort and had an additional 2 months of exposure.  In that light, their cohort efficacy of 79% does not wildly contradict the UK study's lumped efficacy of 88%.
- As far as data going beyond 2.5 months, less than 6% of the UK study would fall into this category compared to roughly 98% of the Israeli study.

Thus, despite the optimism of the UK study, it simply doesn't contain any real information on the efficacies of the vaccine more than 2.5 months after completion, and the data it does show is within 9% (and has overlapping 95% CIs) of the value reported by the Israelis for a similar cohort.  Though the case numbers were small, the Israeli data does contain information on efficacy beyond 4 and 5 months that the UK data cannot currently contradict.

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

<a name="li2021">11</a>: Li B, Deng A, Li K, et al. Viral infection and transmission in a large, well-traced outbreak caused by the SARS-CoV-2 Delta variant. medRxiv. Published online July 23, 2021:2021.07.07.21260122. doi:10.1101/2021.07.07.21260122

<a name="cdcleak">12</a>: Center for Disease Control, [Improving communications around vaccine breakthrough and vaccine effectiveness](https://context-cdn.washingtonpost.com/notes/prod/default/documents/54f57708-a529-4a33-9a44-b66d719070d9/note/7335c3ab-06ee-4121-aaff-a11904e68462.#page=1). Leaked to Washington Post, July 29, 2021

<a name="efficacy_footnote">13</a>: The efficacy cannot be used as a simple probability for determining the odds of a single infection being prevented, since the same individuals being tested more than one time will statistically result in a population incidence ratio less than the efficacy. Nor can it be thought of as a steady-state incidence ratio in this case, since the infection incidence ratio by definition will start as an undefined value 0/0 and then change as time progresses and immunity develops/decays.  

<a name="lopezbernal2021">14</a>: Lopez Bernal J, Andrews N, Gower C, et al. Effectiveness of Covid-19 Vaccines against the B.1.617.2 (Delta) Variant. N Engl J Med. Published online July 21, 2021:NEJMoa2108891. doi:10.1056/NEJMoa2108891

<a name="hall2022">15</a>: Hall VJ, Foulkes S, Charlett A, et al. SARS-CoV-2 infection rates of antibody-positive compared with antibody-negative health-care workers in England: a large, multicentre, prospective cohort study (SIREN). The Lancet. 2021;397(10283):1459-1469. doi:10.1016/S0140-6736(21)00675-9

<a name="israelimoh2021">16</a>: https://www.gov.il/BlobFolder/reports/vaccine-efficacy-safety-follow-up-committee/he/files_publications_corona_two-dose-vaccination-data.pdf 


---

## Model parameter discovery

*This section is a placeholder for descriptions of how model parameters were/are actually estimated from historical data*

### Relative infectivity of delta variant

### Efficacy of vaccines and natural immunity

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
