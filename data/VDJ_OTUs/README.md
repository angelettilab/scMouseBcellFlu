# B-cell receptor repertoire sequencing data

This directory should contain a folder for each sample (e.g., `lung0_2/`, `lung7_1`, `mln7_2`, consistent with the sample naming in the `metadata.csv` file in the parent `data/` directory), with each folder containing the B-cell receptor sequencing annotation file from Cellranger (`filtered_contig_annotations_[SAMPLE].csv`) and the associated FASTA file (`filtered_contig_[SAMPLE].fasta`).

The directory contents should be organized as follows:

```
VDJ_OTUs
├── README.md
├── lung0_2
│   ├── filtered_contig_annotations_lung0_2.csv
│   └── filtered_contig_lung0_2.fasta
├── lung14_2
│   ├── filtered_contig_annotations_lung14_2.csv
│   └── filtered_contig_lung14_2.fasta
├── lung14_3
│   ├── filtered_contig_annotations_lung14_3.csv
│   └── filtered_contig_lung14_3.fasta
├── lung28_1
│   ├── filtered_contig_annotations_lung28_1.csv
│   └── filtered_contig_lung28_1.fasta
├── lung28_2
│   ├── filtered_contig_annotations_lung28_2.csv
│   └── filtered_contig_lung28_2.fasta
├── lung28_3
│   ├── filtered_contig_annotations_lung28_3.csv
│   └── filtered_contig_lung28_3.fasta
├── lung7_1
│   ├── filtered_contig_annotations_lung7_1.csv
│   └── filtered_contig_lung7_1.fasta
├── mln14_1
│   ├── filtered_contig_annotations_mln14_1.csv
│   └── filtered_contig_mln14_1.fasta
├── mln14_2
│   ├── filtered_contig_annotations_mln14_2.csv
│   └── filtered_contig_mln14_2.fasta
├── mln14_3
│   ├── filtered_contig_annotations_mln14_3.csv
│   └── filtered_contig_mln14_3.fasta
├── mln28_1
│   ├── filtered_contig_annotations_mln28_1.csv
│   └── filtered_contig_mln28_1.fasta
├── mln28_2
│   ├── filtered_contig_annotations_mln28_2.csv
│   └── filtered_contig_mln28_2.fasta
├── mln28_3
│   ├── filtered_contig_annotations_mln28_3.csv
│   └── filtered_contig_mln28_3.fasta
├── mln7_1
│   ├── filtered_contig_annotations_mln7_1.csv
│   └── filtered_contig_mln7_1.fasta
├── mln7_2
│   ├── filtered_contig_annotations_mln7_2.csv
│   └── filtered_contig_mln7_2.fasta
├── spleen0_1
│   ├── filtered_contig_annotations_spleen0_1.csv
│   └── filtered_contig_spleen0_1.fasta
├── spleen0_2
│   ├── filtered_contig_annotations_spleen0_2.csv
│   └── filtered_contig_spleen0_2.fasta
├── spleen14_1
│   ├── filtered_contig_annotations_spleen14_1.csv
│   └── filtered_contig_spleen14_1.fasta
├── spleen14_2
│   ├── filtered_contig_annotations_spleen14_2.csv
│   └── filtered_contig_spleen14_2.fasta
├── spleen14_3
│   ├── filtered_contig_annotations_spleen14_3.csv
│   └── filtered_contig_spleen14_3.fasta
├── spleen28_1
│   ├── filtered_contig_annotations_spleen28_1.csv
│   └── filtered_contig_spleen28_1.fasta
├── spleen28_2
│   ├── filtered_contig_annotations_spleen28_2.csv
│   └── filtered_contig_spleen28_2.fasta
├── spleen7_1
│   ├── filtered_contig_annotations_spleen7_1.csv
│   └── filtered_contig_spleen7_1.fasta
└── spleen7_2
    ├── filtered_contig_annotations_spleen7_2.csv
    └── filtered_contig_spleen7_2.fasta
    
```
