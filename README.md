# The macroevolutionary singularity of snakes

This repository is associated with the following article:

Title P.O., Singhal S., Grundler M.C., Costa G.C., Pyron R.A., Colston T.J., Grundler M.R., Prates I., Stepanova N., Jones M.E.H., Cavalcanti L.B.Q., Colli G.R., Di-Poï N., Donnellan S.C., Moritz C., Mesquita D.O., Pianka E.R., Smith S.A., Vitt L.J., Rabosky D.L. 2024. The macroevolutionary singularity of snakes. *Science* 383, 918-923. [link](https://doi.org/10.1126/science.adh2449)

**Looking for the main time-calibrated squamate phylogeny?** Here are shortcut download links for the [main tree](https://raw.githubusercontent.com/macroevolution/squamata/main/data/1.tree_inference/mainTrees/squamates_Title_Science2024_ultrametric_constrained.tre), and the [pseudo-posterior set of trees](https://raw.githubusercontent.com/macroevolution/squamata/main/data/1.tree_inference/pseudoposterior/pseudoposterior.100.trees) (right-click and 'save as').

We have provided coding scripts and data used to infer the squamate phylogeny, prepare trait and other datasets, run comparative analyses, and generate all main and supplemental figures. README files with additional information are provided within repository subfolders.

Please check out our [Dryad repository](https://doi.org/10.5061/dryad.p5hqbzkvb) for a static version of these files that dates to the date of publication.

Although this repository generally reflects the Dryad repository, we will update files here as needed. 

For instance, you will find an **updated version** of the Supplmentary Online Materials, where we have corrected errors and typos. This file is called *SOM-Science-updated.pdf*.

## Description of the data and file structure

The data folder and scripts folder are intended to be placed in a common directory such that file paths in the scripts can point to the appropriate data files.  

Data are organized into several subfolders that are aligned with the different stages of analysis. Scripts generally mirror this structure. 

Code for comparative analyses and main or supplemental figures should be runnable by the user. Code for other stages of analysis are shared for the sake of transparency, but are not runnable as provided. 

**Phylogenies** can be found in the data repository, in data/1.tree_inference > mainTrees

**Diversification and trait data** can be found in data/alldat.csv (see associated readme).

## Changelog

- Update to SOM: Updated/corrected Table S5, which had some incorrect dates (22 Feb 2024)