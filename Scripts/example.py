#!/usr/local/bin/python

from MI_Pack import Atoms
from MI_Pack import Ions
from MI_Pack import Isotopes
from MI_Pack import Identification

if __name__ == "__main__":
    
    ion_mode = "POS"
    atoms = Atoms() #Create standard set of atoms
    ions = Ions(ion_mode) #Create standard set of "positive" ion forms (e.g. H+, Na+, K(39)+)
    isotopes = Isotopes(ion_mode)
    fn_MIDB = 'MIDB_v0_9b.sqlite'

    ID = Identification("peaklist.txt", "peaklist.sqlite", ion_mode, 1.0, atoms, ions, isotopes, fn_MIDB) # Create dataset for identification - 0.0 ppm error
    ID.MIDB("MIDB", {"KEGG_COMPOUND":["hsa"]})
    ID.SPS() # Single-peak search
    ID.EFS() # Empirical formula(e) search
    ID.TM() # Transformation Mapping
    ID.Maps("SPS")
    ID.output("SPS")

                
