#!/usr/bin/env python3
import sys
def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)

#### Import some standard modules
import os
import argparse
import os.path
import timeit
import re
import json
import numpy
import gzip
import pandas

#### Import technical modules and pyteomics
from pyteomics import pepxml, auxiliary


####################################################################################################
#### PepXml Reader class
class PepXmlReader:


    ####################################################################################################
    #### Constructor
    def __init__(self, pepxml_file, verbose=None):

        self.pepxml_file = pepxml_file

        #### Create a place to store our composite spectra for analysis
        self.psms = []

        #### Set verbosity
        if verbose is None: verbose = 0
        self.verbose = verbose


    ####################################################################################################
    #### Read psms
    def read_psms(self):

        #### Set up information
        t0 = timeit.default_timer()
        stats = { 'n_psms': 0 }

        #### Show information
        if self.verbose >= 1:
            eprint(f"INFO: Reading pepXML file {self.pepxml_file}")
            progress_intro = False

        #### If the pepXML is gzipped, then open with zlib, else a plain open
        match = re.search('\.gz$',self.pepxml_file)
        if match:
            infile = gzip.open(self.pepxml_file)
        else:
            infile = open(self.pepxml_file, 'rb')

        #### Define the output columns
        column_name_list = [
            'scan_number',                  # Scan number of the spectrum
            'probability',                  # Probability that the spectrum has been correctly identified
        ]

        rows = []
        psm_number = 0

        #### Read PSMs from the file
        with pepxml.read(infile) as reader:
            for psm in reader:

                #### Set default values
                peptideprophet_probability = None
                iprophet_probability = None
                msrun_name = '?'
                keep = True

                #### Read the values from this PSM
                sequence = psm['search_hit'][0]['peptide']
                charge = psm['assumed_charge']
                spectrum_name = psm['spectrum']

                #### Extract the ms run name from the spectrum name
                match = re.match(r"(.+)\.(\d+)\.(\d+)\.(\d+)$",spectrum_name)
                if match:
                    msrun_name = match.group(1)
                    scan_number = match.group(2)
                    if int(charge) != int(match.group(4)):
                        print(f"ERROR: Charge in spectrum name {spectrum_name} does not match charge {charge}")
                        exit()

                else:
                    print(f"ERROR: Failed to parse spectrum name {spectrum_name}")
                    exit()

                #### Loop over the analysis results, extracting PeptideProphet, iProphet, and PTMProphet probabilities
                mod_probabilities = None
                for analysis_result in psm['search_hit'][0]['analysis_result']:
                    if analysis_result['analysis'] == 'peptideprophet':
                        peptideprophet_probability = analysis_result['peptideprophet_result']['probability']
                    if analysis_result['analysis'] == 'interprophet':
                        iprophet_probability = analysis_result['interprophet_result']['probability']

                #### If this is a PSM that meets the thresholds, then build the USI, and add it to the data rows
                if keep:

                    psm_number += 1

                    row = [ scan_number, peptideprophet_probability ]
                    rows.append(row)


                #### Testing. Print the data structure of one spectrum and exit or just exit early
                #if psm_number ==281:
                #    auxiliary.print_tree(analysis_result['ptmprophet_result'])
                #    sys.exit(10)
                #    auxiliary.print_tree(psm)
                #    sys.exit(10)
                #if stats['n_psms'] > 2000:
                #    break

                #### Update counters and print progress
                stats['n_psms'] += 1
                if self.verbose >= 1:
                    if stats['n_psms']/1000 == int(stats['n_psms']/1000):
                        if not progress_intro:
                            eprint("INFO: Reading psms.. ", end='')
                            progress_intro = True
                        eprint(f"{stats['n_psms']}.. ", end='', flush=True)

        infile.close()
        if self.verbose >= 1: eprint("")

        #### Print final timing information
        t1 = timeit.default_timer()
        eprint(f"INFO: Read {stats['n_psms']} PSMs from {self.pepxml_file}")
        eprint(f"INFO: Elapsed time: {t1-t0}. Processed {stats['n_psms']/(t1-t0)} PSMs per second")

        # Create a pandas data frame from the list data
        table = pandas.DataFrame(rows, columns = column_name_list)

        #### Sort the table by PeptideProphet probability
        #table = table.sort_values(by=['final_probability'], ascending=False)

        #### Write the final dataframe to a TSV file
        output_filename = f"{self.pepxml_file}.tsv"
        table.to_csv(output_filename, sep='\t', index=False)


####################################################################################################
#### For command-line usage
def main():

    argparser = argparse.ArgumentParser(description='Exports a PepXML file to a SiteTab file')
    argparser.add_argument('--min_psm_probability', action='store', default=0.90, help='Minimum PSM (PeptideProphet or iProphet) probability to accept')
    argparser.add_argument('--verbose', action='count', help='If set, print more information about ongoing processing' )
    argparser.add_argument('--version', action='version', version='%(prog)s 0.5')
    argparser.add_argument('filename', type=str, help='Filename of the pepXML file to read')
    params = argparser.parse_args()

    #### Set verbose
    verbose = params.verbose
    if verbose is None:
        verbose = 1

    #### Ensure that the file really is there before starting work
    if not os.path.isfile(params.filename):
        print(f"ERROR: File '{params.filename}' not found or not a file")
        return

    reader = PepXmlReader(params.filename, verbose=verbose)

    reader.read_psms()


#### For command line usage
if __name__ == "__main__": main()
