Small angle scattering programs
===============================

## Overview

The programs in this package were written to process data for small angle x-ray scattering (SAXS) and small angle neutron scattering (SANS).

* calculate_scattering_from_pdb/calc_scat_with_without_params.py (python2)
* calculate_scattering_from_pdb/calc_scat_from_previous_calc.py (python2)
* mulch2_calculate_componenents_from_neutron_contrast_scattering/MULCH2 (C++)
* example_data_neutron_scattering_of_bamlet/BAMLET_neutron_scattering/SANS_and_SAXS_data

#### calculate_scattering_from_pdb/calc_scat_with_without_params.py (python2)

```
This program reads in the scattering curves for each scattering component and cross-term between scattering components, 
and reads in the contrasts and volumes, calculates the new scattering. 
A previous program calculated the scattering curves, and it may have taken a while for that to be calculated. 
This program allows the contrast to be varied and scattering to be recalculated with the varied contrast 
without having to rerun the timeconsuming calculating for component scattering. 
```

#### calculate_scattering_from_pdb/calc_scat_from_previous_calc.py (python2)

```
This python program started of as a python version of a Steen Hanson MCsimul program (written in FORTRAN) 
that was download from http://www.matfys.kvl.dk/~steen/MCsimul.zip May 2012. 
It was changed to read input in the PDB file format. 
It was extended to process multiple densities according to Whitten and Trewhella 2008. 
The Guinier processing in this python program has not yet been tested to ensure that it is still correct after the extensions made to the python program 

This program reads in a file in pdb format, containing points that represent scattering points. 
For x-ray scattering, the points represent electrons. 
This program reads the co-ordinates of the points, 
not the atom type or any other information in the pdb format file that might be used for display. 
The input is not a pdb file describing atomic co-ordinates, 
but can conveniently be used by pdb programs to view the scattering points as if they were atoms. 
The pdb format does not have a field for contrast (or scattering length density). 
This field will be read from columns 81-88. If the field is not filled in the input pdb file, 
then it will be defaulted to the value of 1 for that scattering point. 
```

#### mulch2_calculate_componenents_from_neutron_contrast_scattering/MULCH2 (C++)

```
MULCH2 is a modified version of MULCH software (Whitten et al. 2008) that calculates X-ray and neutron scattering contrasts, radius of gyration from scattering data, and extracts form factors from scattering data. MULCH extracts 2 components. MULCH2 can extract for more than 2 components.
```

#### example_data_neutron_scattering_of_bamlet/BAMLET_neutron_scattering/SANS_and_SAXS_data

```
This data is X-ray and neutron scattering data of the BAMLET compound. MULCH2 and calculate_scattering_from_pdb/calc_scat_with_without_params.py were used on this data in Rath et al. 2017. 
Further processing was carried out as detailed in Rath et al. 2017, and the resulting PDB models are in example_data_neutron_scattering_of_bamlet/BAMLET_neutron_scattering/BAMLET_models
```

## Citations

The scattering programs in this package were published in the following peer-reviewed journal:

Neutron scattering shows a droplet of oleic acid at the center of the BAMLET complex.
Rath EM, Duff AP, Gilbert EP, Doherty G, Knott RB, Church WB.
Proteins. 2017 Jul;85(7):1371-1378. doi: 10.1002/prot.25298. Epub 2017 Apr 22.
PMID: [28380660](https://www.ncbi.nlm.nih.gov/pubmed/28380660) DOI: [10.1002/prot.25298](https://onlinelibrary.wiley.com/doi/abs/10.1002/prot.25298)

MULCh: modules for the analysis of small-angle neutron contrast variation data from biomolecular assemblies
Whitten AE, Caib S, Trewhella J.
Journal of Applied Crystallography. 2008; 41, 222–226.

