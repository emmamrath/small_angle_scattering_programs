

This program is a modified version of the MULCH suite of programs from http://smb-research.smb.usyd.edu.au/NCVWeb/
It has been modified to allow more than 2 components. This was achieved by implementing a singular value decomposition.

To compile the programs :

cp Makefile2 Makefile
rm *.o
make Contrast2
make Rg2
make Compost2

To run the programs on the example 2-component system :

./Contrast2 < contrast2_KinA.txt > contrast2_KinA.html
./Rg2 < rg2_KinA.txt > rg2_KinA.html
./Compost2 < compost2_KinA.txt ./ compost > compost2_KinA.html

The output files are :

contrast2_KinA.html
rg2_KinA.html
compost2_KinA.html
i11_xscat_compost.dat
i22_xscat_compost.dat
i12_xscat_compost.dat
SVD_compost.dat
YdivSigma_compost.dat

To run Compost on a contrived example 3-component system :

./Compost2 < compost2_KinA_3component.txt ./ compost > compost2_KinA_3component.html


