# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pyspedas
import pytplot
import math
import scipy
import numpy
import csv
import os
#import matplotlib


def main():
    pyspedas.erg.pwe_ofa(trange=['2017-03-27', '2017-03-28'], level='l2')
    # pytplot.tplot('erg_pwe_ofa_l2_spec_B_spectra_132')
    freq_data = pytplot.get_data('erg_pwe_ofa_l2_spec_B_spectra_132')
    timestep = freq_data.times
    B_freq = 1000 * freq_data.v  # since the frequency is in kHz
    Pf = freq_data.y
    # Timestep determines each time step of the data from the day (in seconds?)
    # B_freq specifies the frequency of the B-field
    # Pf specifies the P(f) function for the given timestep, which correlates directly to frequency
    #matplotlib.pyplot.plot(timestep, Pf)
    # Check to make sure that the data is of the right shape
    if(len(timestep) != len(Pf)):
        return False

    if(len(B_freq) != len(Pf[0])):
        return False

    # Assuming we're of the right shape, we keep going.
    # We want range [feq/40,feq/2] so to find fce we need to find our B field
    # The following process should find our fce for each 8 second interval for the given day
    # which we then use to find whether the range is acceptable at a minimum
    pyspedas.erg.mgf(trange=['2017-03-27', '2017-03-28'])
    mag_data = pytplot.get_data('erg_mgf_l2_mag_8sec_sm')
    B_field = mag_data.y
    fce = []
    for i in B_field:
        B_mag = math.sqrt(((i[0])**2)+((i[1])**2)+((i[2])**2))
        fce_temp = 29.2*B_mag
        fce.append(fce_temp)

    """    
    03082024
    The above has finished calculating the B_field, which lasts for 8 seconds.
    That 8-second window needs to be checked to see if the condition of fce domain applies
    so we need one loop to cycle through each term in the fce matrix
    and then for each term in the fce matrix, we check the next 8i + a indices,
    where i is the index of the fce matrix
    and a is a constant that increments from 0 to 7 to scroll through the seconds
    of the 8-second window, of the B_freq at the time 8i + a
    to see if that time satisfies.
    If it does, we keep;
    a) the relevant timestep (from which we can calulate the i index of the fce matrix
                              using integer division t = 8i + a // 8)
    b) the B^2 for that timestep
    c) the f_expected
    d) the f_delta
    In a large array. First column is timestep, second column is B^2, third column is f_expected, fourth column is f_delta
    """

    # initializes the output array
    whistler_wave_fingerprints = numpy.array([[0.0, 0.0, 0.0, 0.0]])
    detection_count = 0
    print("start!")
    for i in range(len(fce)):
        for a in range(8):
            index = 8*i + a
            # check for end of sequence
            if index >= 86387:
                break
            [divided_B_freq, divided_Pf_at_t] = subdivide_freq_domain(
                fce[i], B_freq, Pf[index])
            if (isinstance(divided_B_freq, (bool)) or len(divided_Pf_at_t) == 0):
                continue
            elif check_conditions(divided_Pf_at_t) == True:
                [B_freq_prepared, Pf_prepared_at_t] = calc_limits(
                    fce[i], divided_B_freq, B_freq, Pf[index])
                if (len(B_freq_prepared)!=0 and len(Pf_prepared_at_t)!=0):
                    Bw = calc_Bw(Pf_prepared_at_t, B_freq_prepared)
                    braket_f = calc_braket_f(Pf_prepared_at_t, B_freq_prepared, Bw)
                    delta_f = calc_delta_f(Pf_prepared_at_t, B_freq_prepared, braket_f, Bw)
                    if detection_count == 0:
                        whistler_wave_fingerprints[0][0] = timestep[index]
                        whistler_wave_fingerprints[0][1] = Bw
                        whistler_wave_fingerprints[0][2] = braket_f
                        whistler_wave_fingerprints[0][3] = delta_f
                    else:
                        temp = numpy.array([[0.0, 0.0, 0.0, 0.0]])
                        temp[0][0] = timestep[index]
                        temp[0][1] = Bw
                        temp[0][2] = braket_f
                        temp[0][3] = delta_f
                        whistler_wave_fingerprints = numpy.concatenate((whistler_wave_fingerprints, temp), axis=0)
                    detection_count += 1
            else:
                continue
    with open('output.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Time (s)', 'B^2', '<f>', 'del_f'])
        writer.writerows(whistler_wave_fingerprints)
    print("Done!")

    """
    currently, fce is a list, Pf is an array of float32,
    timestep is an array of float64, Bfreq is an array of float32,
    since fce_range is something that'll need to be calculated for each
    of our functions' fce in order to calculate fmin and fstar, then
    we can have a specific function that does just that
    01222024
    
    A function for determining the range of frequencies applicable for our whistler waves
    given the fce of a specific timestep as well as the B_freq of that timestep.
    Does not assume any additional information about the time-step, so the inputs must
    correlate to the correct timestep outside of this function.
    Returns an array of floats containing the frequencies in B_freq which satisfies
    the domain [fce/40, fce/2]. Does not yet include a step to check if there's a full
    domain here. fce is a float, B_freq is a float array
    """


def subdivide_freq_domain(fce, B_freq, Pf):
    upper = fce/2
    lower = fce/40
    B_freq_lower = list(filter(lambda i: i >= lower, B_freq))
    if len(B_freq_lower) != 0:
        B_freq_low = B_freq_lower[0]
    else:
        return[False, False]
    B_freq_higher = list(filter(lambda i: i <= upper, B_freq))
    if len(B_freq_higher) != 0:
        B_freq_high = B_freq_higher[len(B_freq_higher)-1]
    else:
        return[False, False]
    index_low = numpy.where(B_freq == B_freq_low)[0][0]
    index_high = numpy.where(B_freq == B_freq_high)[0][0]
    divided_B_freq = B_freq[index_low:index_high]
    divided_Pf = Pf[index_low:index_high]
    return [divided_B_freq, divided_Pf]


def check_conditions(divided_Pf):
    Pf_max = max(divided_Pf)
    Pf_min = min(divided_Pf)
    if(Pf_max > 3*Pf_min):
        return True
    else:
        return False

# first divides frequency into proper domain, second divides Pf into proper range
# third checks our conditions to keep timestep

# Returns our integral-ready frequencies and densities as our first and second
# indices respectively


def calc_limits(fce, divided_B_freq, B_freq, Pf):
    f_min = min(divided_B_freq)
    f_max = max(divided_B_freq)
    f_star = min([(fce/2), (2*((f_max)-(f_min)))])
    max_entry = list(filter(lambda i: i > f_star, B_freq))
    if len(max_entry) != 0:
        index_high = numpy.where(B_freq == max_entry[0])[0][0]
    else:
        index_high = numpy.where(B_freq == f_max)[0][0]
    min_entry = list(filter(lambda i: i > f_min, B_freq))[
        0]  # What? Double Check 03082024
    index_low = numpy.where(B_freq == min_entry)[0][0]

    """
    so problem; f_star > f_max. Should we just integrate over all of f_max?
    current implementation will include that consideration
    and just integrate over f_max, if needed
    03082024
    """
    Pf_prepared = Pf[index_low:index_high]
    B_freq_prepared = B_freq[index_low:index_high]
    return [B_freq_prepared, Pf_prepared]

# fmin and fstar are floats, Pf is a float array


def calc_Bw(Pf_prepared, B_freq_prepared):
    Bw = scipy.integrate.simpson(Pf_prepared, B_freq_prepared)
    return Bw

# same as above, but with Bw


def calc_braket_f(Pf_prepared, B_freq_prepared, Bw):
    integrand = numpy.multiply(Pf_prepared, B_freq_prepared)
    integral = scipy.integrate.simpson(integrand, B_freq_prepared)
    braket_f = integral/Bw
    return braket_f

# same as above, but with Bw and braket_f
# note that numpy.multiply multipilies the frequency and density elementwise
# and we simpson integrate over the accepted ranges of our frequency and density


def calc_delta_f(Pf_prepared, B_freq_prepared, braket_f, Bw):
    f_factor = []
    for i in B_freq_prepared:
        factor_mag = i*i - braket_f*braket_f
        f_factor.append(factor_mag)
    integrand = numpy.multiply(Pf_prepared, f_factor)
    integral = scipy.integrate.simpson(integrand, B_freq_prepared)
    delta_f = integral/Bw
    return delta_f


main()
