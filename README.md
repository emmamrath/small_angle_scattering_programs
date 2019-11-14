Small angle scattering programs
===============================

## Overview

The programs in this package were written to process data for small angle x-ray scattering (SAXS) and small angle neutron scattering (SANS).

* calculate_scattering_from_pdb/calc_scat_from_previous_calc.py (python2)
* calculate_scattering_from_pdb/calc_scat_with_without_params.py (python2)

#### calculate_scattering_from_pdb

* calc_scat_from_previous_calc.py (python2)

 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; This python program started of as a python version of a Steen Hanson MCsimul program (written in FORTRAN) 
 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; that was download from http://www.matfys.kvl.dk/~steen/MCsimul.zip May 2012. 
 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; It was changed to read input in the PDB file format. 
 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; It was extended to process multiple densities according to Whitten and Trewhella 2008. 
 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The Guinier processing in this python program has not yet been tested to ensure that it is still correct after the extensions made to the python program

 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; This program reads in a file in pdb format, containing points that represent scattering points. 
 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; For x-ray scattering, the points represent electrons. 
 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; This program reads the co-ordinates of the points, 
 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; not the atom type or any other information in the pdb format file that might be used for display. 
 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The input is not a pdb file describing atomic co-ordinates, 
 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; but can conveniently be used by pdb programs to view the scattering points as if they were atoms. 
 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The pdb format does not have a field for contrast (or scattering length density). 
 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; This field will be read from columns 81-88. If the field is not filled in the input pdb file, 
 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; then it will be defaulted to the value of 1 for that scattering point. 

* calc_scat_with_without_params.py (python2)

 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; This program reads in the scattering curves for each scattering component and cross-term between scattering components, 
 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; and reads in the contrasts and volumes, calculates the new scattering. 
 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; A previous program calculated the scattering curves, and it may have taken a while for that to be calculated. 
 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; This program allows the contrast to be varied and scattering to be recalculated with the varied contrast 
 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; without having to rerun the timeconsuming calculating for component scattering. 

## Citation

The scattering programs in this package were published in the following peer-reviewed journal:

Neutron scattering shows a droplet of oleic acid at the center of the BAMLET complex.
Rath EM, Duff AP, Gilbert EP, Doherty G, Knott RB, Church WB.
Proteins. 2017 Jul;85(7):1371-1378. doi: 10.1002/prot.25298. Epub 2017 Apr 22.
PMID: [28380660](https://www.ncbi.nlm.nih.gov/pubmed/28380660) DOI: [10.1002/prot.25298](https://onlinelibrary.wiley.com/doi/abs/10.1002/prot.25298)
