#!/bin/python2

import wave_util
import os
import sys

def zerodown_data(argv):
    # if a data directory is specified from the command line, use it
    if len(argv) > 1:
        datpath = argv[1]
    # if not, look in the directory the program was launched from
    else:
        curpath = os.path.dirname(os.path.abspath(__file__))
        datpath = os.path.join(curpath, 'data')

    # if a sampling rate is specified use it, if not use a default value of 0.5
    if len(argv) > 2:
        sr = argv[2]
    else:
        sr = 0.5

    # lists of relevant hourly data
    hmean = []
    hmax = []
    hhts = []
    hht10 = []

    # loop through all our data files and compute relevant hourly data
    """
    The weakness to this approach is that individual waves spanning files are not identified (they could but it would be a hassle). It also means we need to handle border cases (zero-dimension matrices) in our zerodown function.
    """
    for file in os.listdir(datpath):
        #if file == '0065':
        # read raw data from file
        with open(os.path.join(datpath, file), 'r') as f:
            raw = f.read().split('\n')

        # remove first line (header) and any empty lines
        raw = [r for r in raw[1:] if r]
        # retrieve elevation data
        data = [r.split(',')[1].replace('-\.','-0\.') for r in raw]
        # zerodowncross
        mean, max, hts, ht10 = wave_util.zerodown(data, sr)
        # add to our hourly list
        #print file, mean, max, hts, ht10
        hmean.append(mean)
        hmax.append(max)
        hhts.append(hts)
        hht10.append(ht10)

    # write results to file
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'mean_data'), 'w') as f:
        [f.write(str(mean[0]) + ',' + str(mean[1]) + '\n') for mean in hmean]
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'max_data'), 'w') as f:
        [f.write(str(max[0]) + ',' + str(max[1]) + '\n') for max in hmax]
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'hts_data'), 'w') as f:
        [f.write(str(hts[0]) + ',' + str(hts[1]) + '\n') for hts in hhts]
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ht10_data'), 'w') as f:
        [f.write(str(ht10[0]) + ',' + str(ht10[1]) + '\n') for ht10 in hht10]

if __name__ == ('__main__'):
    zerodown_data(sys.argv)

