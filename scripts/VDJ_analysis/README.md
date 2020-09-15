# BCR Seq Analysis

VDJ (B-cell receptor) sequencing data was processed using the [Immcantation toolbox](https://immcantation.readthedocs.io/en/stable/), and later merged with the RNA-Seq data in the Seurat object (or AnnData object in Scanpy).

## VDJ workflow

The entire VDJ processing workflow can be executed using the `immcantation_VDJ_workflow.sh` bash script within this directory. Be sure to first specify the `main` project directory (`main=...`) at the top of the script.

```
cd scripts/VDJ_analysis/
bash immcantation_VDJ_workflow.sh
```

The following sections describe each step of the VDJ workflow included in the `immcantation_VDJ_workflow.sh` bash script.


### Define directories

Change the `main` directory such that it corresponds to the path on your local system. The script uses this path to find and create other directories and data relative to `main`, and assumes that the directory structure matches that of this GitHub repository.

```
# specify main project directory
main='/MY/PATH/scMouseBcellFlu'
```

### Prepare environment

Activate the **immcant-env** conda environment containing the packages necessary to run the pipeline. If the `immcant-env` conda environment has not been created, the conda environment `.yml` files are in the [`envs/`](../../envs) directory of this repository, and can be used to create the required conda environments (more details can be found within the [`envs/`](../../envs) directory).

```
source activate immcant-env
```

The Immcantation pipeline requires some accessory scripts that are not included in a conda or python package, but instead must be retrieved from [their Bitbucket repository](https://bitbucket.org/kleinstein/immcantation/src/master/scripts/). The necessary files have already been retrieved from their repository and placed in the [`scripts/immcantation/`](../immcantation) subdirectory, so it is only necessary to add this subdirectory to your PATH so bash is able to find the scripts within.

```
# On MacOS, for example
PATH="$PATH:"$main"/scripts/immcantation"
```


### Retrieve Ig BLAST reference databases

The Immcantation accessory functions are used to retrieve and build the mouse Ig reference sequences. The output of these commands is placed in the [`data/immcantation/`](../../data/immcantation) subdirectory.

```
fetch_igblastdb.sh -o $main/data/immcantation/igblast
fetch_imgtdb.sh -o $main/data/immcantation/germlines/imgt
imgt2igblast.sh -i $main/data/immcantation/germlines/imgt -o $main/data/immcantation/igblast
```

*Note that the imgt2igblast.sh script uses the `realpath` command that does not come natively installed on MacOS, and therefore may need to be installed using e.g. homebrew:*

```
brew install coreutils
```

### Run Ig BLAST on VDJ FASTA files

Use the `AssignGenes.py` script to BLAST the VDJ FASTA files and assign reads to Ig genes. The `*` wildcards act to process fasta files for all samples within the `data/VDJ_OTUs/` subdirectory.


### Parse VDJ files into Change-O database format

Process and filter each VDJ seq sample individually to produce Change-O database (`.tab`) files. These `.tab` files are tab-delimited files containing information on each read, such as V/D/J assignments, their actual sequences, etc. This parsing will also split each sample file into heavy (IGH) and light (IGK, IGL) chain files.


### Merge samples by mouse and estimate sequence distance threshold

Run the `01_VDJ_group_mice_and_threshold.R` script to group sequence data by mouse and calculate the distribution of Hamming distances between IGH V segment sequences for each mouse. The script will automatically predict the threshold Hamming distance that defines which sequences should be treated as distinct clones (those with a distance above the threshold), and which should be assigned to the same clone (those with a distance below the threshold).

The results of the Hamming distance distributions and threshold estimation are placed in mouse-specific subdirectories (`M1`, `M2`, `M3`, etc.) within the `analysis/immcantation/clone_assignment/` directory, where predicted thresholds are drawn as a vertical dashed red line, and summarized in the enclosed `predicted_thresholds.csv` file.

Note that for some mice, the automatic threshold prediction fails (because the distribution is not bimodal) and therefore a default value (0.1) is assumed. This default threshold can be changed using the `--default_threshold` flag input to the `01_VDJ_group_mice_and_threshold.R` script. To manually assign thresholds based on another method or visual inspection, simply edit the `predicted_thresholds.csv` file.


### Define clones and germline sequences

With the VDJ sequence Hamming distance thresholds estimated in the previous step, heavy chains are assigned a clone ID using the `DefineClones.py` script.

The `light_cluster.py` script is used to update the resulting clone assignments based on light chain(s) that are paired to each heavy chain (inferred based on matching CELL identifiers in each file). This script also removes any cells that are associated with more than one heavy chain.

The `CreateGermlines.py` script then infers the germline sequence of the V segment corresponding to each heavy and light chain, and appends this to the sequence database files.


### Quantify VDJ mutation burden

The R-script `02_VDJ_mutation_quant.R` is used to estimate the mutation burden of each heavy and light chain sequence based on its difference from its inferred germline sequence.


### Merge heavy and light chain data and add to Seurat object

All of the VDJ sequence, clonotype, and mutation data obtained from the steps above can be integrated with the scRNA-seq Seurat object (in the `metadata` slot) using the `03_VDJ_RNAseq_integration.R` script. Heavy and light chains are paired using the CELL barcodes, where multiple light chains can be associated with each heavy chain.

The resulting Seurat object allows for integrated analysis of VDJ and mRNA sequencing data; for example by overlaying VDJ mutation frequencies onto the transcriptome-based UMAP.







