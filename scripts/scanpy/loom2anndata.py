import os, sys, argparse
import scanpy as sc


parser = argparse.ArgumentParser(description='Convert loom to scanpy anndata')
parser.add_argument('-i', '--input', type=str, help='Input loom file')
parser.add_argument('-o', '--output', type=str, help='Output h5ad file')

args = parser.parse_args()


#in_loom = "/Users/asbj/projects/sc_projects/single-cell-hackathon-2020/datasets/bone_marrow/scanpy/10x/filt_seurat_object.loom"
#out_adata = "/Users/asbj/projects/sc_projects/single-cell-hackathon-2020/datasets/bone_marrow/scanpy/10x/filt_seurat_object.h5ad"

adata = sc.read_loom(args.input)
adata.write(args.output)


