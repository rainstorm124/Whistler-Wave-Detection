# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 16:18:53 2024

@author: raini
"""

import scipy, numpy

def main():
    y = numpy.arange(0,10)
    x = [0.1, 0.3, 0.5, 0.7, 0.9, 1.1]
    value = scipy.integrate.simpson(y, x)
    return value

main()