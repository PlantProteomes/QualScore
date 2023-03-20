#!/usr/bin/env python3

import sys
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

#### Import some standard modules
import os
import argparse
import os.path
import timeit
import re
import json
import numpy
import pickle
import gzip
from pyteomics import mzml, auxiliary


####################################################################################################
#### mzML Assessor class
class MzMLReader:


    ####################################################################################################
    #### Constructor
    def __init__(self, mzml_file, verbose=None):
        self.mzml_file = mzml_file

        #### Create a place to store our spectra for analysis
        self.spectra = {}

        #### Set verbosity
        if verbose is None: verbose = 0
        self.verbose = verbose


    ####################################################################################################
    #### Read spectra
    def read_spectra(self, write_fragmentation_type_file=None):

        #### Set up information
        t0 = timeit.default_timer()
        stats = {
            'n_spectra': 0,
            'n_ms0_spectra': 0,
            'n_ms1_spectra': 0,
            'n_ms2_spectra': 0,
        }

        #### Show information
        if self.verbose >= 1:
            eprint(f"\nINFO: Assessing mzML file {self.mzml_file}", flush=True)
            progress_intro = False

        #### If the mzML is gzipped, then open with zlib, else a plain open
        if self.mzml_file.endswith('.gz'):
            infile = gzip.open(self.mzml_file)
        else:
            infile = open(self.mzml_file, 'rb')

        #### Read spectra from the file
        with mzml.read(infile) as reader:
            #try:
            if True:
                for spectrum in reader:

                    #### Debugging. Print the data structure of the first spectrum
                    #if stats['n_spectra'] == 0:
                    #    auxiliary.print_tree(spectrum)
                    #    print(spectrum)
                    #    return

                    #### Extract the MS level of the spectrum
                    ms_level = spectrum['ms level']

                    #### If the ms level is 2, then examine it for information
                    if ms_level == 2 and 'm/z array' in spectrum:
                        precursor_mz = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']
                        #### Try to get the charge information
                        try:
                            charge_state = int(spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state'])
                        except:
                            charge_state = 'unknown'

                        peaklist = {
                            'm/z array': spectrum['m/z array'],
                            'intensity array': spectrum['intensity array']
                        }

                        #### Check for zero length and very sparse spectra
                        if len(spectrum['m/z array']) == 0:
                            # zero length spectrum!
                            pass

                    #### Update counters and print progress
                    stats['n_spectra'] += 1
                    if stats['n_spectra']/500 == int(stats['n_spectra']/500):
                       eprint(f"{stats['n_spectra']}.. ", end='', flush=True)
            #except:
            else:
                self.log_event('ERROR','MzMLCorrupt',f"Pyteomics threw an error reading mzML file! File may be corrupt. Check file '{self.mzml_file}'")

        infile.close()
        #if self.verbose >= 1: eprint("")

        #### Print final timing information
        t1 = timeit.default_timer()
        print(f"\nINFO: Read {stats['n_spectra']} spectra from {self.mzml_file} in {int(t1-t0)} sec ({stats['n_spectra']/(t1-t0)} spectra per sec)", end='', flush=True)



####################################################################################################
#### For command-line usage
def main():

    argparser = argparse.ArgumentParser(description='Read an mzML file')
    argparser.add_argument('--verbose', action='count', help='If set, print more information about ongoing processing' )
    argparser.add_argument('files', type=str, nargs='+', help='Filenames of one or more mzML files to read')
    params = argparser.parse_args()

    #### Set verbose
    verbose = params.verbose
    if verbose is None: verbose = 1

    #### Loop over all the files to ensure that they are really there before starting work
    for file in params.files:
        if not os.path.isfile(file):
            print(f"ERROR: File '{file}' not found or not a file")
            return

    #### Loop over all the files processing them
    for file in params.files:

        #### Assess the mzML file
        reader = MzMLReader(file, verbose=verbose)
        reader.read_spectra()


#### For command line usage
if __name__ == "__main__": main()
