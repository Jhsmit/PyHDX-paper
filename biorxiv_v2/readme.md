README
======

The scripts in this directory can be used to reproduce the figures in the bioRxiv manuscript
by Smit et al: 'Probing universal protein dynamics using residue-level Gibbs free energy' 
(DOI: https://doi.org/10.1101/2020.09.30.320887 )


Installation
============

Create and activate the conda environment

```console
$ conda env create -f environment.yml 
$ conda activate p38_biorxiv_v2 
```

Download and install PyHDX  

```
$ pip install PyHDX==0.4.0b1
``` 

Data availability
=================

The datasets for the proteins ecSecB and SecA (https://doi.org/10.1016/j.str.2021.03.015) are included in this repository.

Datasets published elsewhere:  

mtSecB: https://doi.org/10.1038/s41467-019-08747-4
hPREP:  https://doi.org/10.1038/s41598-017-02550-1 

The datasets for PpiA, PpiB, TF and MBP will be part of upcoming publications.
This means that currently the code will not run out-of-the-box due to missing files. We apologize for the inconvenience.

Running
=======

Run the numbered python scripts in order to run PyHDX analysis. 
Next run 0x_run_all_figures to generate all figures.