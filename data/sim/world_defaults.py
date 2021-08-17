"""
    This module contains functionality to build default properties for the world based on the current literature and
    parameter estimation
"""
from typing import List

import numpy
import scipy.stats
from scipy import interpolate
from sim.dynamics import WorldProperties, VariantProperties, DiscreteFunction


def default_world_properties() -> WorldProperties:
    alpha = VariantProperties(
        incubation=alpha_variant_incubation(),
        infectivity=alpha_infectivity_curve(),
        vax_immunity=pfizer_alpha_efficacy(),
        natural_immunity=natural_alpha_efficacy()
    )

    delta = VariantProperties(
        incubation=delta_variant_incubation(),
        infectivity=delta_infectivity_curve(),
        vax_immunity=pfizer_delta_efficacy(),
        natural_immunity=natural_delta_efficacy()
    )

    return WorldProperties(alpha, delta)


def alpha_variant_incubation() -> List[float]:
    """
    The alpha/pre-alpha variant incubation curve comes from "Early Transmission Dynamics in Wuhan, China, of Novel
    Coronavirus–Infected Pneumonia" by Li et al.

    The value is returned as a set of probabilities extracted from the cumulative density function of the exponential
    weibull used in Li et al. The paper does not contain the parameters, so they had to be discovered by curve fitting
    to values extracted from an image capture of Figure 2.

    The return value is a set of cumulative probabilities in a list where the list index is the day and the value is
    the cumulative probability of that day or earlier having been selected.
    :return: a list of cumulative probabilities for each day's index
    """
    # These values came from a least squares minimization of points taken off of Li et al.
    a, c, scale = [12.6245844, 0.57604443, 0.64524305]

    days = []
    for i in range(21):
        days.append(scipy.stats.exponweib.cdf(i, a, c, 0, scale))

    return days


def delta_variant_incubation() -> List[float]:
    """
    The estimation for the delta variant's incubation curve is simply compressing the alpha variant down to 2/3rds of
    its domain, based on the 4 day to 6 day ratio observed in "Viral infection and transmission in a large, well-traced
    outbreak caused by the SARS-CoV-2 Delta variant".

    :return: a list of cumulative probabilities for each day's index
    """
    # These values came from a least squares minimization of points taken off of Li et al.
    a, c, scale = [12.6245844, 0.57604443, 0.64524305]

    days = []
    for i in range(21):
        scaled = i * 3 / 2.0  # Scale the domain variable larger to shrink the space
        days.append(scipy.stats.exponweib.cdf(scaled, a, c, 0, scale))

    return days


def alpha_infectivity_curve() -> DiscreteFunction:
    """
    The alpha variant's infectivity curve is estimated in Ashcroft et al. (doi:10.4414/smw.2020.20336) by a gamma
    distribution with the properties below.
    """
    # Ashcroft et al., this distribution is more centered around the day of symptom onset
    shape, rate, shift = 97.18750, 3.71875, 25.6250

    # The original He et al. curve, which leans more towards pre-symptomatic transmission. Ashcroft claimed to have
    # found a mistake in the original He paper, but I did not look for He's response.
    # shape, rate, shift = 1.5625, 0.53125, 2.12500

    # Generate an inverse cdf table
    offset = 10
    total_days = 20
    values = []
    for i in range(total_days):
        y = scipy.stats.gamma.pdf(i - offset, shape, -shift, 1 / rate)
        values.append(y)

    # The last value in the infectivity curve MUST be 0, otherwise the simulation will not be able to efficiently
    # remove infected individuals after they stop being contagious
    values.append(0)

    return DiscreteFunction(offset, values)


def custom_infectivity_curve(shape: float, rate: float, shift: float) -> DiscreteFunction:
    """
    Produces a custom infectivity curve based on a gamma distribution
    """
    limit = 0.0005

    # Start by walking backwards until we hit the limit
    values = []
    median = round(scipy.stats.gamma.median(shape, -shift, 1 / rate))
    position = median
    traveled = 0

    value = scipy.stats.gamma.pdf(position, shape, -shift, 1 / rate)
    while value > limit:
        values.append(value)
        traveled += 1
        position -= 1
        value = scipy.stats.gamma.pdf(position, shape, -shift, 1 / rate)

    values.reverse()

    position = median + 1
    value = scipy.stats.gamma.pdf(position, shape, -shift, 1 / rate)
    while value > limit:
        values.append(value)
        position += 1
        value = scipy.stats.gamma.pdf(position, shape, -shift, 1 / rate)

    values[0] = 0
    values[-1] = 0
    return DiscreteFunction(traveled - median - 1, values)

    # Generate an inverse cdf table
    offset = 10
    total_days = 20
    values = []
    for i in range(total_days):
        y = scipy.stats.gamma.pdf(i - offset, shape, -shift, 1 / rate)
        values.append(y)

    # The last value in the infectivity curve MUST be 0, otherwise the simulation will not be able to efficiently
    # remove infected individuals after they stop being contagious
    values.append(0)

    return DiscreteFunction(offset, values)


def delta_infectivity_curve() -> DiscreteFunction:
    """
    The delta variant's infectivity curve is currently estimated by scaling the alpha curve by a factor of 3
    """
    return alpha_infectivity_curve().scale_y(2.0)


def natural_alpha_efficacy() -> DiscreteFunction:
    """
    The best data on reinfection I could find was from a study in the UK that followed a healthcare provider cohort
    which got tested early.  From that I could at least make a rough estimate of reinfection.

    SARS-CoV-2 infection rates of antibody-positive compared with antibody-negative health-care workers in
    England: a large, multicentre, prospective cohort study (SIREN).

    doi:10.1016/S0140-6736(21)00675-9

    :return:
    """
    # The number of cases between the previously seropositive cohort and the negative cohort at the end of the
    #
    positive_cohort = 153 / 8278
    negative_cohort = 1704 / 17383

    # The median interval of reinfection for individuals with date data for the first infection was 201 days, but
    # ranged between 95 and 257 days. The study included data from June 18, 2020 to Jan 11, 2021. For reasonably
    # estimable cases the range spanned from 90 days to 345 with a median at 241.
    #
    # Suffice it to say that reinfection within 90 days is improbable, but after that immunity wanes to some degree
    # allowing cases to accumulate until they reach the odds ratio observed.

    adjust_day = 200
    efficacy = (negative_cohort - positive_cohort) / negative_cohort
    estimated_ratio = [(0, 1), (90, 1), (adjust_day, efficacy), (300, efficacy)]
    nat = interpolate.interp1d(*zip(*estimated_ratio))
    values = []
    for i in range(300):
        values.append(float(nat(i)))
    return DiscreteFunction(0, values).mean_filter(50)


def natural_delta_efficacy() -> DiscreteFunction:
    """
    From "SARS-CoV-2 variants of concern and variants under investigation in England, Technical Briefing 19" dated
    July 23, 2021:

        The adjusted odds ratio of reinfection with the Delta variant was 1.46 (95% CI 1.03 to
        2.05) compared to the Alpha variant. The risk of reinfection was not elevated for Delta if
        the primary infection was <180 days (adjusted odds ratio = 0.79, 95% CI 0.49 to1.28) but
        was higher for those with a prior infection ≥180 days earlier (adjusted odds ratio = 2.37,
        95%CI 1.43 to 3.93). Further work to examine the risk of reinfection is being undertaken.

    Based on this, I'm going to do the same thing as the alpha variant except the target efficacy will be adjusted so
    that the adjusted odds ratio matches the 2.37 number

    Todo: there should be a delta->delta efficacy, as it should matter which variant you get
    :return:
    """
    # The positive and negative cohorts from the past study
    positive_cohort = 153 / 8278
    negative_cohort = 1704 / 17383

    adjust_day = 200
    efficacy = (negative_cohort - (positive_cohort * 2.37)) / negative_cohort
    estimated_ratio = [(0, 1), (150, 1), (adjust_day, efficacy), (300, efficacy)]
    nat = interpolate.interp1d(*zip(*estimated_ratio))
    values = []
    for i in range(300):
        values.append(float(nat(i)))

    return DiscreteFunction(0, values).mean_filter(50)


def pfizer_delta_efficacy() -> DiscreteFunction:
    """

    :return:
    """

    return pfizer_alpha_efficacy().scale_y(0.9)


def pfizer_alpha_efficacy() -> DiscreteFunction:
    """
    Based on "Six Month Safety and Efficacy of the BNT162b2 MRNA COVID-19 Vaccine", (doi:10.1101/2021.07.28.21261159),
    which is still in pre-print and not peer reviewed.
    """

    # Because the paper doesn't have their raw data, I had to use https://apps.automeris.io/wpd/ to extract it

    vaccinated = [(0, 0),
                  (5.978021978021985, 0.0010052356020942288),
                  (7.659340659340671, 0.0013916230366492118),
                  (8.967032967032978, 0.0016753926701570665),
                  (20.17582417582419, 0.0020104712041884715),
                  (46.32967032967032, 0.002178010471204181),
                  (60.90109890109889, 0.0025968586387434406),
                  (89.48351648351648, 0.003350785340314133),
                  (109.47252747252746, 0.004020942408376957),
                  (127.78021978021975, 0.005361256544502618),
                  (153.56043956043953, 0.007036649214659685),
                  (174.85714285714283, 0.007706806282722509),
                  (190.5494505494505, 0.009047120418848156)]

    placebo = [(0, 0),
               (4.109890109890117, 0.0006701570680628238),
               (29.890109890109898, 0.007036649214659685),
               (48.94505494505495, 0.01156020942408377),
               (77.34065934065934, 0.02043979057591623),
               (89.29670329670328, 0.02580104712041885),
               (109.47252747252746, 0.034848167539267005),
               (123.67032967032968, 0.041047120418848164),
               (130.021978021978, 0.04523560209424082),
               (153.56043956043953, 0.054785340314136115),
               (176.72527472527474, 0.06299476439790574),
               (188.30769230769226, 0.06567539267015705),
               (194.28571428571428, 0.0691937172774869)]

    vax_x, vax_y = zip(*vaccinated)
    vax = interpolate.interp1d(vax_x, vax_y, kind="linear")

    p_x, p_y = zip(*placebo)
    pl = interpolate.interp1d(p_x, p_y, kind="linear")

    last_day = int(min(max(vax_x), max(p_x)))
    x_ = numpy.linspace(0, last_day)
    efx = []
    for i in range(last_day + 1):
        if pl(i) <= 0:
            efx.append(0)
        else:
            efx.append(float(max(0, (pl(i) - vax(i)) / pl(i))))

    # import matplotlib.pyplot as plt
    # plt.plot(x_, efx)
    # plt.show()

    return DiscreteFunction(0, efx).mean_filter(10)
