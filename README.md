# POLARe 

Differential elastic scattering cross-section of spin-polarized positron and electron beams impinging on neutral atoms

This code allows to calculate the differential elastic scattering cross-section (DESCS) of positrons and electrons impinging on neutral atoms.
It uses the relativistic partial wave expansion method (Mott cross-section).
The z axis is along the direction of particle incidence.

The user should include the input data in the "INPUT/input.txt" file. For example, for calculating the DESCS of 1000 eV electrons in Au (atomic number Z=79), the "INPUT/input.txt" file can contain the following information:

Z 79
E(eV) 1000
Positron/Electron(0/1) 1 
Px 0
Py 0
Pz 1
Azimuth(deg) 0

As the z axis was chosen along the direction of particle incidence, with these input data the spin-polarization of the beam is longitudinal. When the beam is spin-polarized with transverse polarization, i.e. when Px and/or Py are set to be different from zero, the DESCS also depends on the azimuthal angle (and on the spin-polarization as well). 

The "INPUT/screening.txt" file contaions all the parameters necessary for calculating the screening function. From Z=1 to Z=18: Best fit of the Hartree-Fock atomic model (H.L. Cox Jr. and R.A. Bonham, The Journal of Chemical Physics, Vol. 47, pp. 2599-2608, 1967). From Z=19 to Z=54: Best fit of the Dirac-Hartree-Fock-Slater atomic model (H.L. Cox Jr. and R.A. Bonham, The Journal of Chemical Physics, Vol. 47, pp. 2599-2608, 1967). From Z=55 to Z=92: Best fit of the Dirac-Hartree-Fock-Slater atomic model (F. Salvat, J.D. Martinez, R. Mayol, J. Parellada, Physical Review A, Vol. 36, pp. 467-474, 1987). 

For electron beams, the exchange effects are included, according to Salvat and Mayol (F. Salvat and R. Mayol, Computer Physics Communications, Vol. 74, pp. 358-374, 1993), using the Furness and McCarthy semiphenomenological optical model (J.B. Furness and I.E. McCarthy, Journal of Physics B: Atomic and Molecular Physics, Vol. 6, pp. 2280- 229 1, 1973).

The results are saved in the folder RESULTS and are graphically represented in the folder FIGURES.
