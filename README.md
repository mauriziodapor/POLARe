# POLARe 

# Differential elastic scattering cross-section of spin-polarized positron and electron beams impinging on neutral atoms

POLARe is a C code that allows you to calculate the differential elastic scattering cross section (DESCS) of positrons and electrons hitting neutral atoms [M. Dapor, Journal of Physics B: Atomic, Molecular and Optical Physics 55, 095202 (2022)].

It uses the relativistic method of partial wave expansion (Mott theory).

Since the z-axis was chosen along the direction of particle incidence, the spin polarization of the beam at Pz = 1 is longitudinal. If the beam is spin-polarized with transverse polarization, i.e., if Px and/or Py are different from zero, the DESCS depends on both the azimuthal angle and the spin polarization.

The user should record the input data in the file "INPUT/input.txt". For example, to calculate the DESCS of 1000 eV electrons in U (atomic number Z = 92) with transverse spin polarization Px = 1 and azimuthal angle = 90 degrees, the "INPUT/input.txt" file should contain the following information:

------------------------

Z 92

E(eV) 1000

Positron/Electron(0/1) 1

Px 1

Py 0

Pz 0

Azimuth(deg) 90

------------------------

The "INPUT/screening.txt" file contains all the parameters required to calculate the screening function. 

From Z=1 to Z=17: Best fit of the Hartree-Fock atomic model [H.L. Cox Jr. and R.A. Bonham, The Journal of Chemical Physics 47, 2599 (1967)]. 

From Z=18 to Z=54: Best fit of the Dirac-Hartree-Fock-Slater atomic model [H.L. Cox Jr. and R.A. Bonham, The Journal of Chemical Physics 47, 2599 (1967)]. 

From Z=55 to Z=92: Best fit of the Dirac-Hartree-Fock-Slater atomic model [F. Salvat, J.D. Martinez, R. Mayol, J. Parellada, Physical Review A 36, 467 (1987)]. 

In order to be able to use the screening function of Salvat et al. also in the range of atomic numbers Z = 1-54, you must rename the file "INPUT/Salvatetal1987.txt" to "INPUT/screening.txt".

For electron beams, the exchange effects according to Salvat and Mayol [F. Salvat and R. Mayol, Computer Physics Communications 74, 358 (1993)] are considered using the semiphenomenological optical model of Furness and McCarthy [J.B. Furness and I.E. McCarthy, Journal of Physics B: Atomic and Molecular Physics 6, 2280 (1973)].

The results (DESCS for both spin-polarized and spin-unpolarized beams, S, T and U parameters, atomic density, and cumulative probability) are stored in the RESULTS folder and displayed graphically in the FIGURES folder.

Comparisons of selected POLARe calculations of the DESCS of spin-unpolarized electron beams (atomic numbers: 18, 29, 36, and 79; electron energies: 1keV and 4keV) with the calculations of Riley et al. [M.E. Riley, J. MacCallum, and F. Biggs, Atomic Data and Nuclear Data Tables 15, 443 (1975)] can be found in the COMPARISONS folder. 

Details on POLARe (theory and numerical approach) can be found in the book:

Maurizio Dapor. Electron–Atom Collisions: Quantum-Relativistic Theory and Exercises, Berlin, Boston: De Gruyter, 2022. https://doi.org/10.1515/9783110675375

Compile with the following command:

gcc -O2 -g POLARe.c -o POLARe -lgsl -lgslcblas -lm

Execute it with the following command:

./POLARe

# POLARe2

# Differential elastic scattering cross-section of spin-polarized positron and electron beams impinging on neutral atoms (with correlation polarization potential)

Note that POLARe2 contains the correlation polarization potential according to F. Salvat, Physical Review A 68, 012708 (2003).

Long-range polarization for electrons and positrons: Buckingham potential.

Short-range polarization potential (electrons): J.P. Perdew and A. Zunger, Physical Review B 23, 5048 (1981).

Short-range polarization potential (positrons): A. Jain, Physical Review A 41, 2437 (1990).

Compile with the following command:

gcc -O2 -g POLARe2.c -o POLARe2 -lgsl -lgslcblas -lm

Execute it with the following command:

./POLARe2

# SPAS

# Spin-Polarization After Scattering

The program SPAS enables the calculation of the spin polarization of electron and positron beams after elastic scattering [M. Dapor, Physics Open 14, 100134 (2023)]. 

It requires as input the functions S, T and U (which can be obtained with POLARe or POLARe2), the initial spin polarization, and the azimuthal angle.

Compile with the following command:

gcc -O2 -g SPAS.c -o SPAS -lgsl -lgslcblas -lm

Execute it with the following command:

./SPAS
