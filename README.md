# POLARe 

Differential elastic scattering cross-section of spin-polarized positron and electron beams impinging on neutral atoms

POLARe is a c code which allows to calculate the differential elastic scattering cross-section (DESCS) of positrons and electrons impinging on neutral atoms.

It uses the relativistic partial wave expansion method (Mott theory).

The user should include the input data in the "INPUT/input.txt" file. For example, for calculating the DESCS of 1000 eV electrons in U (atomic number Z=92), the "INPUT/input.txt" file can contain the following information:

Z 92

E(eV) 1000

Positron/Electron(0/1) 1

Px 1

Py 0

Pz 0

Azimuth(deg) 90

As the z axis was chosen along the direction of particle incidence, when Pz=1 the spin-polarization of the beam is longitudinal. When the beam is spin-polarized with transverse polarization, i.e. when Px and/or Py are set to be different from zero, the DESCS depends on the azimuthal angle and on the spin-polarization as well. 

The "INPUT/screening.txt" file contains all the parameters necessary for calculating the screening function. From Z=1 to Z=17: Best fit of the Hartree-Fock atomic model [H.L. Cox Jr. and R.A. Bonham, The Journal of Chemical Physics 47, 2599 (1967)]. From Z=18 to Z=54: Best fit of the Dirac-Hartree-Fock-Slater atomic model [H.L. Cox Jr. and R.A. Bonham, The Journal of Chemical Physics 47, 2599 (1967)]. From Z=55 to Z=92: Best fit of the Dirac-Hartree-Fock-Slater atomic model [F. Salvat, J.D. Martinez, R. Mayol, J. Parellada, Physical Review A 36, 467 (1987)]. In order to use the screening function of Salvat et al. even in the range of atomic numbers Z = 1-54, it is necessary to rename the file "INPUT/Salvatetal1987.txt" which must become "INPUT/screening.txt".

For electron beams, the exchange effects are included, according to Salvat and Mayol [F. Salvat and R. Mayol, Computer Physics Communications 74, 358 (1993)], using the Furness and McCarthy semiphenomenological optical model [J.B. Furness and I.E. McCarthy, Journal of Physics B: Atomic and Molecular Physics 6, 2280 (1973)].

The results (DESCS for both spin-polarized and spin-unpolarized beams, S, T, and U parameters, atomic density, and cumulative probability) are saved in the folder RESULTS and are graphically represented in the folder FIGURES.

Comparisons of POLARe selected calculations of the DESCS of spin-unpolarized electron beams (atomic numbers: 18, 29, 36, and 79; electron energies: 1keV and 4keV) with the calculations of Riley et al. [M.E. Riley, J. MacCallum, and F. Biggs, Atomic Data and Nuclear Data Tables 15, 443 (1975)] are presented in the folder COMPARISONS. 

Details about POLARe (theory and numerical approach) can be found in the book:

Maurizio Dapor. Electron–Atom Collisions: Quantum-Relativistic Theory and Exercises, Berlin, Boston: De Gruyter, 2022. https://doi.org/10.1515/9783110675375

Compile using the following command:

gcc -O2 -g POLARe.c -o POLARe -lgsl -lgslcblas -lm

Run using the following command:

./POLARe

Note that POLARe2 includes the correlation-polarization potential according to F. Salvat, Physical Review A 68, 012708 (2003).
Long-range polarization for both electrons and positrons: Buckingham potential.
Short-range polarization potential (electrons): J.P. Perdew and A. Zunger, Physical Review B 23, 5048 (1981).
Short-range polarization potential (positrons): A. Jain, Physical Review A 41, 2437 (1990).

Compile using the following command:

gcc -O2 -g POLARe2.c -o POLARe2 -lgsl -lgslcblas -lm

Run using the following command:

./POLARe2

The SPAS program (Spin-Polarization After Scattering) allows to calculate the spin-polarization of electron and positron beams after elastic scattering. It requires in input the functions S, T, and U (which can be obtained using POLARe or POLARe2), the initial spin-polarization, and the azimuthal angle.

Compile using the following command:

gcc -O2 -g SPAS.c -o SPAS -lgsl -lgslcblas -lm

Run using the following command:

./SPAS
