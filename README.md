# pyFermi
![Fermi surface of Aluminum, color encoded by its spin-related properties <b>2](newplot.png)

## Introduction
This is an web application that could reconstruct the fermi surface and some related physical properties from ab initio calculation, such as VASP (Vienna ab-initio simulation packages), quantum espresso.

The reconstruction process use Wannier function to interpolate the fermi surface and wavefunctions. Formalism from this paper:

Yates, J. R., Wang, X., Vanderbilt, D., & Souza, I. (2007). Spectral and Fermi surface properties from Wannier interpolation. Physical Review B - Condensed Matter and Materials Physics, 75(19), 1–12. https://doi.org/10.1103/PhysRevB.75.195121

The visualization and interaction are achieved by using DASH.

Currently in its early stage.

