#!/usr/bin/env python
# released under the GNU General Public License version 3.0 (GPLv3)

from distutils.core import setup
import glob 
import sys
import os

def main():
    if sys.version_info[0] != 2 and sys.version_info[1] != 7:
        raise Error("Python-2.7.8 is required")
    example_files = glob.glob(os.path.join('Scripts', '*.*'))
    example_path = os.path.join('Scripts', 'MI_Pack')

    setup( name="MI-Pack",
    version= "2.0 (Galaxy)",
    description= "Metabolite Identification Package (MI-Pack)",
    author= "Ralf Weber, Mark Viant",
    author_email= "weberrj@bham.ac.uk, m.viant@bham.ac.uk",
    url= "http://www.biosciences.bham.ac.uk/labs/viant/",
    license= "GPL",
    platforms = ['Windows XP, UNIX'],
    keywords = ['Metabolomics', 'Mass spectrometry', 'metabolites', 'identification', 'Data mining', 'Mass errors'],
    packages=['MI_Pack'],

    package_dir={'MI-Pack': 'MI_Pack'},
    package_data={'MI_Pack': ['*.pyd', '*.so']},
     
    data_files=[(example_path, example_files)],
    long_description="MI-Pack - Metabolite Identification Package",)

if __name__ == "__main__":
    main()

