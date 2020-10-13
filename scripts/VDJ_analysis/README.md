# BCR Seq Analysis

VDJ (B-cell receptor) sequencing data was processed using the [Immcantation toolbox](https://immcantation.readthedocs.io/en/stable/), and later merged with the RNA-Seq data in the Seurat object (or AnnData object in Scanpy).

## VDJ workflow

The entire VDJ processing workflow can be executed using the `immcantation_VDJ_workflow.sh` bash script within this directory. Be sure to first specify the main project directory (`main=...`) at the top of the script.

```
cd scripts/VDJ_analysis/
bash immcantation_VDJ_workflow.sh
```

Note that there is a `metadata.csv` file in the [`data/`](../../data) directory of the repository which contains information about each of the samples to be processed. For example:

| dataset | assay | chemistry | mouse_nr | infection | day_post_infection | organ | organ_day |
| ------- | ----- | --------- | -------- | --------- | ------------------ | ----- |---------- |
| spleen0_1 | rna | v2 | M7 | naive | D0 | spleen | spleen0 |
| spleen0_2 | rna | v2 | M8 | naive | D0 | spleen | spleen0 |
| spleen7_1 | rna | v2 | M9 | infected | D7 | spleen | spleen7 |
| spleen7_2 | rna | v2 | M10 | infected | D7 | spleen | spleen7 |
| spleen14_1 | rna | v2 | M1 | infected | D14 | spleen | spleen14 |
| spleen14_2 | rna | v2 | M2 | infected | D14 | spleen | spleen14 |
| spleen14_3 | rna | v2 | M3 | infected | D14 | spleen | spleen14 |
| ... | ... | ... | ... | ... | ... | ... | ... |

If a different dataset is to be used, this `metadata.csv` file will need to be updated accordingly. It is important that the **`dataset`** column is present and the entries match the names of the folders in the [`data/VDJ_OTUs/`](../../data/VDJ_OTUs) directory. Also, the **`mouse_nr`** column must be present in order to combine datasets based on subject (mouse) when performing the clone assignment.

The following sections describe each step of the VDJ workflow included in the `immcantation_VDJ_workflow.sh` bash script.


### Define directories

Change the main directory such that it corresponds to the path on your local system

```
# specify main project directory
main='/Users/your/path/here/d_angeletti_1910'
```

### Prepare environment

Activate the **immcant-env** conda environment containing the packages necessary to run the pipeline. If the `immcant-env` conda environment has not been created, the conda environment `.yml` files are in the [`envs/`](../../envs) directory of this repository, and can be used to create the required conda environments (more details can be found within the [`envs/`](../../envs) directory).

```
source activate immcant-env
```

The Immcantation pipeline requires some accessory scripts that are not included in a nicely packaged conda or python package, but instead must be retrieved from [their Bitbucket repository](https://bitbucket.org/kleinstein/immcantation/src/master/scripts/). The necessary files have already been retrieved from their repository and placed in the [`scripts/immcantation/`](../immcantation) subdirectory, so it is only necessary to add this subdirectory to your PATH so bash is able to find the scripts within.

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

Process and filter each VDJ seq sample individually to produce Change-O database (`.tab`) files. These `.tab` files are just tab-delimited files containing information on each read, such as V/D/J assignments, their actual sequences, etc.


### Infer Ig genotype and estimate sequence distance threshold

Run the `01_VDJ_genotype_and_threshold.R` script to retrieve (unmutated) germline VDJ sequences and to calculate the distribution of Hamming distances between measured sequences for each mouse. The script will automatically predict the threshold Hamming distance that defines which mutants should be treated as distinct clones (those with a distance above the threshold), and which should be collapsed into the same germline sequence (those with a distance below the threshold).

The results of the Hamming distance distributions and threshold estimation are placed in the `analysis/immcantation/threshold_estimation/` subdirectory, where predicted thresholds are drawn as a vertical dashed red line, and summarized in the `predicted_thresholds.csv` file.

*Note that for some mice, the automatic threshold prediction fails (because the distribution is not bimodal) and therefore a default value (0.1) is assumed.*


### Define clones and germline sequences

Using the VDJ sequence Hamming distance thresholds estimated in the previous step, each read is assigned to a germline sequence to define individual clones.


### Quantify VDJ mutation burden

The R-script `02_VDJ_mutation_quant.R` is used to estimate the mutation burden of each sequence based on its difference from its inferred germline sequence.


### Add VDJ data to Seurat object

All of the VDJ sequence, clonotype, and mutation data obtained from the steps above can be integrated with the scRNA-seq Seurat object (in the `metadata` slot) using the `03_VDJ_RNAseq_integration.R` script. This allows for integrated analysis of VDJ and mRNA sequencing data; for example by overlaying VDJ mutation frequencies onto the transcriptome-based UMAP.







