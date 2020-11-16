# scRNA-seq Cellranger output files

This directory should contain a folder for each sample (e.g., `lung0_2/`, `lung7_1`, `mln7_2`, consistent with the sample naming in the `metadata.csv` file in the parent `data/` directory), with each folder containing the following outputs from the cellranger pipeline: `barcodes.tsv`, `genes.tsv`, `matrix.mtx`.

The directory contents should be organized as follows:

```
cellranger
├── README.md
├── lung0_2
│   ├── barcodes.tsv
│   ├── genes.tsv
│   └── matrix.mtx
├── lung14_2
│   ├── barcodes.tsv
│   ├── genes.tsv
│   └── matrix.mtx
├── lung14_3
│   ├── barcodes.tsv
│   ├── genes.tsv
│   └── matrix.mtx
├── lung28_1
│   ├── barcodes.tsv
│   ├── genes.tsv
│   └── matrix.mtx
├── lung28_2
│   ├── barcodes.tsv
│   ├── genes.tsv
│   └── matrix.mtx
├── lung28_3
│   ├── barcodes.tsv
│   ├── genes.tsv
│   └── matrix.mtx
├── lung7_1
│   ├── barcodes.tsv
│   ├── genes.tsv
│   └── matrix.mtx
├── mln14_1
│   ├── barcodes.tsv
│   ├── genes.tsv
│   └── matrix.mtx
├── mln14_2
│   ├── barcodes.tsv
│   ├── genes.tsv
│   └── matrix.mtx
├── mln14_3
│   ├── barcodes.tsv
│   ├── genes.tsv
│   └── matrix.mtx
├── mln28_1
│   ├── barcodes.tsv
│   ├── genes.tsv
│   └── matrix.mtx
├── mln28_2
│   ├── barcodes.tsv
│   ├── genes.tsv
│   └── matrix.mtx
├── mln28_3
│   ├── barcodes.tsv
│   ├── genes.tsv
│   └── matrix.mtx
├── mln7_1
│   ├── barcodes.tsv
│   ├── genes.tsv
│   └── matrix.mtx
├── mln7_2
│   ├── barcodes.tsv
│   ├── genes.tsv
│   └── matrix.mtx
├── spleen0_1
│   ├── barcodes.tsv
│   ├── genes.tsv
│   └── matrix.mtx
├── spleen0_2
│   ├── barcodes.tsv
│   ├── genes.tsv
│   └── matrix.mtx
├── spleen14_1
│   ├── barcodes.tsv
│   ├── genes.tsv
│   └── matrix.mtx
├── spleen14_2
│   ├── barcodes.tsv
│   ├── genes.tsv
│   └── matrix.mtx
├── spleen14_3
│   ├── barcodes.tsv
│   ├── genes.tsv
│   └── matrix.mtx
├── spleen28_1
│   ├── barcodes.tsv
│   ├── genes.tsv
│   └── matrix.mtx
├── spleen28_2
│   ├── barcodes.tsv
│   ├── genes.tsv
│   └── matrix.mtx
├── spleen7_1
│   ├── barcodes.tsv
│   ├── genes.tsv
│   └── matrix.mtx
└── spleen7_2
    ├── barcodes.tsv
    ├── genes.tsv
    └── matrix.mtx

```
