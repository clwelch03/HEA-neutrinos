# HEA-neutrinos

## Description

Python and Mathematica implementations of the opacities found in [arXiv:2601.18051](https://arxiv.org/abs/2601.18051). Computes total and differential opacities for neutrinos and antineutrinos in hot matter below nuclear density. All internal calculations are done in natural units in powers of fm.  For examples of usage and code to generate the plots found in the manuscript, see `llmcplots.ipynb`. The Mathematica implementation additionally performs a numerical calculation of opacities with a Monte Carlo summation over Landau levels. The results of such a calculation are in the folder `./src/opacity csvs`. See also [the repo](https://github.com/clwelch03/NS-landau-quantization) associated with [our prior publication.](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.111.063009)

## Citation

See `CITATION.cff` for an example citation. 

## AI use statement

AI was used in the preparation of this manuscript for fixing syntax in the data tables and for generating BibTeX entries. It was not used to perform any computations or to contribute to the text of the manuscript.The authors do not consent to the use of any text, data, or code associated with this publication as training data for any AI model.

## License

BSD 3-Clause License
See `LICENSE` file for details.