# slab-dl-monte

Fix bug for doing GCMC in slab/slit geometry in DL_MONTE.

## The problem

For GCMC simulation in the slab geometry where the slit is placed in the center of the box and vacuum on the side,
the density of the fluid at the center of the slit does not match with the bulk density at the same chemical potential.


The reason for this is because in GCMC, the volume of the box $V$ enters the acceptance criteria for insertion as

```math
\mathrm{acc}(o\rightarrow n) = \mathrm{min}\left[ 1, \frac{V}{(N+1)\Lambda^3}\exp(-\beta \Delta U - \mu )\right]
```
and deletion as
```math
\mathrm{acc}(o\rightarrow n) = \mathrm{min}\left[ 1, \frac{N \Lambda^3}{V}\exp(-\beta \Delta U + \mu )\right]
```

In the original implemetation, when the slab geometry is used (with `zfrac` enabled), insertion moves are not attempted in 
the vacuum. This has the consequence of effectively reducing the volume of the system, but this was not accounted for in the acceptance criteria.
As a consequence, there is a knock-on effect of (effectively) changing the chemical potential of the system, so the densities
inside the slit did not correspond to their bulk values. 


This is fixed by changing the code to attempt to insert atoms/molecules everywhere in the box, and then reject/accept accordingly.
Essentially, replacing instances of  
```buff(3) = duni()*(slit_zfrac2-slit_zfrac1) + slit_zfrac1```
to 
```buff(3) = duni() - 0.5_wp```
in the `cell_module.f90`.

## Tests

We tested this fix with GCMC simulations of a Lennard Jones fluid confined between two soft walls (inputs are provided in the `tests` folders).

<p align="center">
<img src="https://github.com/user-attachments/assets/ff685834-f464-4aff-bda1-b2819b407fe3" width="700">
</p>



## Usage
Instruction to install and compile DL_MONTE is at [DL_MONTE official repository](https://gitlab.com/dl_monte/DL_MONTE-2).

To fix the bug, simply replace the `cell_module.f90` file in the source code to the `cell_module.f90` file of this repository. Then recompile and run simulations.

Follow the issue raised on the development gitlab: [issue](https://gitlab.com/dl_monte/DL_MONTE-2/-/issues/30).



