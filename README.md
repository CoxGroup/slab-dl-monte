# slab-dl-monte

Fix bug for doing GCMC in slab/slit geometry in DL_MONTE.

## The problem

For GCMC simulation in the slab geometry where the slit is placed in the center of the box and vacuum on the side,
the density of the fluid at the center of the slit does not match with the bulk density at the same chemical potential.

The reason for this is because insertion moves are not attempted in the vacuum, but the acceptance criteria is still given as 

```math
\mathrm{acc}(o\rightarrow n) = \mathrm{min}\left[ 1, \frac{V}{(N+1)\Lambda^3}\exp(-\beta \Delta U - \mu )\right]
```
where $V$ is the volume of the box. Therefore, the simulation does not sample the correct Boltzmann distribution. 

This is fixed by changing the code to attempt to insert atoms/molecules everywhere in the box, and then reject/accept accordingly.
Specifically, replacing 
`buff(3) = duni()*(slit_zfrac2-slit_zfrac1) + slit_zfrac1` 
with
`buff(3) = duni() - 0.5_wp`
in the `cell_module.f90`.

## Usage
Instruction to install and compile DL_MONTE is at [DL_MONTE official repository](https://gitlab.com/dl_monte/DL_MONTE-2)

To fix the bug, simply replace the `cell_module.f90` in the source code to the `cell_module.f90` file of this repository. Then recompile and run simulations.



