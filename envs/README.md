# envs

This directory contains files defining [conda environments](https://docs.conda.io/projects/conda/en/latest/index.html) and [singularity containers](https://sylabs.io/guides/3.5/user-guide/) used for different functions and packages (corresponding to the names of the subdirectories in this folder).


## Conda environments on different operating systems

Note that conda has problems with compatibility between different operating systems (Windows, MacOS, Linux, etc.).

### MacOS

The conda environments in this project were primarily designed for macOS, so they should work as provided. 

### Windows/PC

Many conda packages are unfortunately not available or not regularly updated for Windows. As a result, attempting to create a conda environment from the `.yml` environment file will often fail due to the required packages not being found.

If you have Windows 10, a workaround is to use the built-in Linux virtual machine (VM) to build and use the conda environment. This VM is called the Windows Subsystem for Linux (WSL). Further instructions on using conda with WSL can be found [here](https://nbisweden.github.io/workshop-scRNAseq/conda_instructions.html) (scroll down to the Windows-specific section in that link).

Another option is to instead use a compute cluster, such as UPPMAX. The operating system in such cases will often be Linux. See the Linux-specific section below for further information.

### UNIX/Linux

Although many conda packages are similarly available for MacOS and Linux, differences exist that can cause the environment build to fail. If this occurs, try editing the environment `.yml` files by removing the version number after each package. For example, change
```
 - r-seurat=3.0.2
```
to
```
 - r-seurat
```

**IMPORTANT!** If you are using a compute cluster (such as UPPMAX), file number and size limitations may prevent you from using conda on these systems. This is often due to the large numbers of files required to define a conda environment. If you are running into performance issues or exceeding file/size limits on a cluster, it is recommended that you instead use *Singularity* (see next section).


## Singularity

[Singularity](https://sylabs.io/guides/latest/user-guide/) wraps up code and environment into a single container file (`.sif`). This means that it will function identically regardless of the operating system it is run on. _However_, Singularity currently only works on UNIX/Linux systems, though there is a beta MacOS version with some limited functionality.

Many of the subdirectories in this directory contain pre-built `.sif` files that can be used on e.g., UPPMAX instead of a conda environment. There are also instructions on how to re-build the `.sif` file, but note that this will require that you have a Mac or Linux with Singularity installed.

To run, for example, a bash script within a singularity container on UPPMAX, use the following command:
```
singularity run my/path/file.sif bash my/path/bash_script.sh
```

Where `my/path/` is replaced with the actual path to the `.sif` file and your script file, and `file.sif` is the actual name of the image file (e.g., `sauron.sif`).


## Subdirectories in this folder

### immcantation

- `immcant-environment.yml` is the conda environment definition file used to create a conda environment (named "immcant-env") that contains the R packages in the [Immcantation framework](https://immcantation.readthedocs.io/en/stable/). Specifically, it contains the `alakazam`, `shazam`, and `tigger` packages. Create the environment using the following command:

```
conda env create -f immcant-environment.yml  # creates environment named "immcant-env"
conda activate immcant-env  # activate the environment
```


### sauron

- `sauron.sif` is the singularity image file for the Sauron package, containing all dependencies and packages necessary to run Sauron functions.

- `Singularity_remote.def` is the singularity definition file used to build the `sauron.sif` image file.

- `build_singularity_image.sh` is a short bash script that contains instructions on how to build the `sauron.sif` image file using the `Singularity_remote.def` file.

NOTE: The conda environment file for the Sauron package can be found in the [Sauron GitHub repository](https://github.com/NBISweden/sauron).


### rna_velocity

- `rvelo_environment.yml` is the conda environment definition file used to create a conda environment (named "rvelo-env") that contains the R dependencies required to perform RNA velocity analysis with the R implementation of [velocyto](http://velocyto.org/).

- `scvelo_environment.yml` is the conda environment definition file used to create a conda environment (named "scvelo") that contains the python packages necessary to perform RNA velocity and trajectory analysis using the [scVelo](https://scvelo.readthedocs.io/api.html) python package. 

- `velocyto_environment.yml` is the conda environment definition file used to create a conda environment (named "velocyto-env") that contains the packages required to run the [velocyto](http://velocyto.org/) command-line tool for generating spliced and unspliced `.loom` files from the scRNA-Seq `.bam` files.

- `velocyto.sif` is the singularity image file containing the velocyto command-line environment (velocyto-env), which has all the dependencies and packages necessary to run velocyto.

- `velocyto_singularity.def` is the singularity definition file used to build the `velocyto.sif` image file.


### trajectory

- `trajectory.sif` is the singularity image file for the R [slingshot](https://bioconductor.org/packages/release/bioc/vignettes/slingshot/inst/doc/vignette.html) and [tradeSeq](https://www.bioconductor.org/packages/release/bioc/vignettes/tradeSeq/inst/doc/tradeSeq.html) trajectory inference packages, containing all dependencies and packages necessary to run these tools.

- `trajectory_singularity.def` is the singularity definition file used to build the `trajectory.sif` image file.

- `build_singularity_image.sh` is a short bash script that contains instructions on how to build the `trajectory.sif` image file using the `trajectory_singularity.def` file.

- `trajectory_environment.yml` is the conda environment definition file used to create a conda environment (named "trajectory-env") that contains the dependencies required to use the [slingshot](https://bioconductor.org/packages/release/bioc/vignettes/slingshot/inst/doc/vignette.html) and [tradeSeq](https://www.bioconductor.org/packages/release/bioc/vignettes/tradeSeq/inst/doc/tradeSeq.html) R packages.








