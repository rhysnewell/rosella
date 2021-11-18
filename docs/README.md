
Rosella
=============

![Rosella logo](/images/rosella.png)
![](https://travis-ci.com/rhysnewell/rosella.svg?branch=master)
![](https://anaconda.org/bioconda/rosella/badges/license.svg)
![](https://anaconda.org/bioconda/rosella/badges/version.svg)
![](https://anaconda.org/bioconda/rosella/badges/platforms.svg)

Rosella is a metagenomic binning algorithm using UMAP and HDBSCAN. It is written in Rust with a python component that 
handles calls to UMAP and HDBSCAN. Rosella aims to be as user friendly as possible with multiple usage modes and installation
methods. 

Please note that Rosella is under active development with new commits often providing much improved results. If you would like
the most up to date version of Rosella please pull the code from `dev` branch. Hopefully releases will stabilise very soon.

## Additional resources

Rosella makes use of a couple new and daunting algorithms. UMAP in particular is an amazing algorithm but might be cause 
for concern since it is difficult to understand how it works and what it is doing. So please look over this amazing article 
by Andy Coenen and Adam Pearce: [Understanding UMAP](https://pair-code.github.io/understanding-umap/)

## Citation

*Watch this space* A paper is on its way. If you use rosella and like the results before the paper, then please cite this GitHub

## License

Code is [GPL-3.0](LICENSE)