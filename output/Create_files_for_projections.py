#!/usr/bin/env python
# coding: utf-8

# Create the files needed for projecting external data onto the [Stemformatics integrated atlases](stemformatics.org/atlas/myeloid). These files were used to generate the figures in the Stemformatics paper.
# 
# Two external datasets were used: 
# - Monkley et.al. RNAseq characterisation of iPSC-derived monocytes, macrophages and dendritic cells and their their primary blood derived counterparts (DOI: 10.1371/journal.pone.0243807)
# - Rosa et.al. Single-cell transcriptional profiling informs efficient reprogramming of human somatic cells to cross-presenting dendritic cells (DOI: 10.1126/sciimmunol.abg5539).
# 
# Note that data files are not included in this repository due to size limits, but the providence of each input file used here 
# is described below.

# In[8]:


# Set up the environment.
import pandas, scanpy


# ### Process the Monkley data

# In[3]:


"""Process the Monkley data. The original gene expression matrix and sample metadata were downloaded from ArrayExpress:
https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-9670. 
Gene expression file: iPSCd_comb_geneTPM.txt
Sample metadata file: E-MTAB-9670.sdrf.txt
All other files are processed and produced here.
"""
def processMonkley():
    # Show what original data look like
    df = pandas.read_csv("data/Monkley/iPSCd_comb_geneTPM.txt", sep="\t", index_col=0)
    samples = pandas.read_csv("data/Monkley/E-MTAB-9670.sdrf.txt", sep="\t", index_col=0)
    display(df.head())
    display(samples.head())
    print("Shape of expression matrix:", df.shape, "Shape of sample matrix:", samples.shape)

    # Sample matrix has more samples than there are in the expression matrix.
    # Create files where the samples are matched between the two.
    
    commonSampleIds = set(df.columns).intersection(set(samples.index))
    df = df[list(commonSampleIds)]
    df.to_csv("data/Monkley/iPSCd_comb_geneTPM.SampleMatched.tsv", sep="\t")
    
    # We also don't need all the columns in the sample matrix - select a few we want and rename them too.
    print(samples.columns.tolist())
    samples = samples[['Characteristics[cell_type_abbrev]','Characteristics[cell type]','Characteristics[cell line]','Characteristics[stimulus]']]
    samples = samples.loc[df.columns].fillna('')
    samples = samples[~samples.index.duplicated(keep='first')]
    samples.columns = [item.replace('Characteristics[','').replace(']','') for item in samples.columns]
    
    # Insert a new column in samples called "combined" which combines cell type with stimulus
    samples.insert(0, 'combined', [f"{row['cell_type_abbrev']}_{row['stimulus']}".replace("_none","") for index,row in samples.iterrows()])
    
    print("Shape of expression matrix:", df.shape, "Shape of sample matrix:", samples.shape)
    display(samples.head())
    samples.to_csv("data/Monkley/samples.tsv", sep="\t")
    
processMonkley()


# In[4]:


"""For the purpose of creating figures for the paper, we further filtered out some samples for simplicity, 
as it can be hard to see patterns with too many cell types on a static printed figure (it works much better
as an interactive plot on the website when above files are used directly to project the data).
"""
def createMonkleySubset():
    samples = pandas.read_csv("data/Monkley/samples.tsv", sep="\t", index_col=0)
    keep = ['HMDM_D3_unstim_63_S21', 'HMDM_D2_unstim_60_S18', 'HMDM_D1_unstim_57_S15', 'IPSDDC_1D_d21_P3_9_S7', 
            'IPSDDC_2D_d32_30_S28', 'IPSDDC_2D_d22_28_S26', 'IPSDDC_3D_d32_47_S7', 'IPSDDC_2D_d29_29_S27', 
            'IPSDDC_1D_d35_P4_13_S11', 'IPSDDC_1D_d28_P3_11_S9', 'IPSDDC_3D_d22_45_S5', 'IPSDDC_3D_d29_46_S6', 
            'IPSDM1M_rep1_03_S3', 'IPSDM_3B_d21_unstim_34_S32', 'IPSDM3M_rep2_18_S18', 'IPSDM3M_rep3_19_S19', 
            'IPSDM_2B_d28_unstim_20_S18', 'IPSDM2M_rep3_11_S11', 'IPSDM_3B_d35_unstim_40_S38', 'IPSDM2M_rep1_09_S9', 
            'IPSDM_2B_d21_unstim_17_S15', 'IPSDM1M_rep2_04_S4', 'IPSDM_2B_d35_unstim_23_S21', 
            'IPSDM_3B_d28_unstim_37_S35', 'IPSDM3M_rep1_17_S17', 'IPSDM2M_rep2_10_S10', 'CD14N_1D_d35_4_S3', 
            'CD14N_1D_d31_3_S2', 'CD14N_1D_d28_2_S1', 'CD14N_2D_d43_26_S24', 'CD14P2_07_S7', 'CD14P3_15_S15', 
            'CD14N_3D_d43_43_S3', 'CD14P1_01_S1', 'CD14P_3B_d21_31_S29', 'CD14P_3D_d43_44_S4', 'CD14M2_08_S8', 
            'CD14P_2B_d21_14_S12', 'CD14P_3B_d28_32_S30', 'CD14P_1D_d31_7_S5', 'CD14P_2B_d35_16_S14', 
            'CD14M3_16_S16', 'CD14P_1D_d35_8_S6', 'CD14P_1D_d28_6_S4', 'CD14P_3B_d35_33_S31', 'CD14P_2B_d28_15_S13', 
            'CD14P_2D_d43_27_S25', 'CD14M1_02_S2', 'HMDDC_D2_unstim_53_S11', 'HMDDC_D3_unstim_55_S13', 
            'HMDDC_D1_unstim_51_S9', 'PBMC_D3_68_S26', 'PBMC_D1_66_S24', 'PBMC_D2_67_S25']
    print("Do specified sample ids all belong to the sample matrix?", set(keep).issubset(set(samples.index)))
    # - in case there was a mistake in this step
    
    samples = samples.loc[keep]
    print("Filtered sample matrix shape:", samples.shape)
    
    # Write the expression and sample matrices which match these samples to files.
    # These were the files used to make projections onto the Stemformatics dendritic cell atlas 
    # (stemformatics.org/atlas/dc). Screenshots were then used to create the figures for the paper.
    df = pandas.read_csv("data/Monkley/iPSCd_comb_geneTPM.SampleMatched.tsv", sep="\t", index_col=0)
    df = df[samples.index]
    df.to_csv("data/Monkley/iPSCd_comb_geneTPM.SampleMatched.subset.tsv", sep="\t")
    samples.to_csv("data/Monkley/samples.subset.tsv", sep="\t")
    
createMonkleySubset()


# ### Process the Rosa data
# We started by downloading GSE162650_hef_all_counts.txt from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162650). Then we ran the following R code to convert this to a sparse matrix, matrix.mtx.
# Other files created by this script include genes.tsv (gene symbols corresponding to the features of the matrix),
# barcodes.tsv (barcodes corresponding to the samples of the matrix) and celltypes.tsv (cell type that each sample
# belongs to).
# 
# ```R
# .libPaths('~/libs/R_libs')
# require(data.table)
# require(Matrix)
# dat <- fread('./GSE162650_hef_all_counts.txt')
# barcode <- colnames(dat)[-1]
# celltype <- substr(barcode,1,nchar(barcode)-19)
# unique_celltype <- unique(celltype)
# gene <- as.vector(as.matrix(dat[,1]))
# dat <- dat[,-1]
# dat <- Matrix(as.matrix(data.frame(dat)),sparse = TRUE)
# writeMM(obj = dat, file = paste(".","/matrix.mtx",sep = ''))
# write(x = gene, file = paste(".","/genes.tsv",sep = ''))
# write(x = barcode, file = paste(".","/barcodes.tsv",sep = ''))
# write(x = celltype, file = paste(".","/celltypes.tsv",sep = ''))
# ```
# 
# We then used scanpy's read_mtx function to read matrix.mtx, and assigned the genes from genes.tsv as its features (vars) and barcodes from barcodes.tsv as its samples (obs). We also log normalised using scanpy's log1p function (to base 2). 
# This created LogNormalised.h5ad file.
# 

# In[2]:


"""Function to output a pseudo-bulk aggregated data, given single cell data as input. 
> pseudoBulkAggregate(df, ['HSC','B','B',...])

df (pandas.DataFrame): genes (rows) x samples (columns) expression matrix (single cell RNA-seq dataset)
groups (pandas.Series): specifies what group each sample belongs to (eg. cell type); index equals df.columns
n (int): how many points to represent the smallest cluster with

This function first determines the smallest group based on groups, then tries to split that group into n sub-groups.
This determines the cluster size which will be used to aggregate. So if the smallest group has 30 samples in it,
for example, and n=3, 10 samples are used to create each aggregated sample for each group. Hence keep an eye on
the sample size, as it may not work well if this is too small.
"""
def pseudoBulkAggregate(df, groups, n=3, verbose=True):
    import random
    
    valueCounts = groups.value_counts().sort_values()
    smallest = valueCounts.tolist()[0] # size of the smallest cluster
    if verbose: print("value counts of cells in the group", valueCounts)
    
    if smallest<=3*n: # if the smallest cluster isn't at least 3 times the value of n, just use all of that cluster
        n = 1
    sampleSize = int(smallest/n)  # how many samples to put into each aggregated cluster
    if sampleSize==1:  # no aggregation happening
        raise Exception("Could not determine how many samples to aggregate together. Likely due to a cluster having only one member.")
    if verbose: print("sample size for each cluster:", sampleSize)

    aggregated = {}
    for item in valueCounts.index:  # each item of group
        subset = df[groups[groups==item].index]  # expression for this sampleGroup
        # we need new sample ids based on sample group eg. "Mono1__0"
        i = 0
        columns = subset.columns.tolist()
        while len(columns)>sampleSize:
            selectedColumns = random.sample(list(columns), sampleSize)
            aggregated["%s__%s" % (item, i)] = subset[selectedColumns].sum(axis=1)
            columns = set(columns).difference(set(selectedColumns))
            i += 1

    return pandas.DataFrame(aggregated, index=df.index)


# In[3]:


"""Read the LogNormalised.h5ad file from above and apply the aggregation function to create aggregated data.
"""
def aggregateRosa():
    adata = scanpy.read("data/Rosa/LogNormalised.h5ad")
    celltypes = open("data/Rosa/celltypes.tsv").read().split("\n")
    celltypes = pandas.Series([item for item in celltypes if item!=''], index=adata.obs.index)

    agg = pseudoBulkAggregate(adata.to_df().transpose(), celltypes)
    samples = pandas.DataFrame({'celltype':[item.split("__")[0] for item in agg.columns]}, index=agg.columns)
    print(agg.shape)
    display(agg.head())
    agg.to_csv("data/Rosa/Rosa_aggregated.tsv", sep="\t")
    samples.to_csv("data/Rosa/Rosa_aggregated_samples.tsv", sep="\t")
    
aggregateRosa()


# In[4]:


# Rosa data comes with gene symbols in the expression matrix - change these to gene ids.
# We downloaded gene id to gene symbol mapping from Ensembl Biomart and named the file EnsemblGenes_HomoSapiens_v91.txt.
def useEnsemblGeneIds():
    ensembl = pandas.read_csv("data/Ensembl/EnsemblGenes_HomoSapiens_v91.txt", sep="\t", index_col=0)
    display(ensembl.head())
    geneIdsFromSymbol = {}
    for geneId,row in ensembl.iterrows():
        if row['gene_name'] not in geneIdsFromSymbol: geneIdsFromSymbol[row['gene_name']] = []
        geneIdsFromSymbol[row['gene_name']].append(geneId)
    
    df = pandas.read_csv("data/Rosa/Rosa_aggregated.tsv", sep="\t", index_col=0)
    print(df.shape)
    rows,geneIds = [],[]
    for geneSymbol,row in df.iterrows():
        for geneId in geneIdsFromSymbol.get(geneSymbol, []):
            rows.append(row.tolist())
            geneIds.append(geneId)
    df = pandas.DataFrame(rows, index=geneIds, columns=df.columns)
    print(df.shape)
    df.to_csv("data/Rosa/Rosa_aggregated_EnsemblId.tsv", sep="\t")
    
useEnsemblGeneIds()


# In[6]:


# Aggregated data frame above was still too large for upload to Stemformatics. Reduce by selecting some cell types.
# These were the final files used to project onto the Stemformatics Myeloid Atlas and the results are shown in the paper.
def createSmallerAggregatedData():
    samples = pandas.read_csv("data/Rosa/Rosa_aggregated_samples.tsv", sep="\t", index_col=0)
    print(samples.shape)
    celltypes = ['DC1','DC2','PDC','Day3','Day6','Day9','HEF']
    samples = samples.loc[[item for item in samples.index if item.split("_")[0] in celltypes and item.split("_")[1]=='B'\
                          or item.startswith('Day9_DP_B')]]
    df = pandas.read_csv("data/Rosa/Rosa_aggregated_EnsemblId.tsv", sep="\t", index_col=0)
    df = df[samples.index]
    print(samples.shape, df.shape)
    df.to_csv("data/Rosa/Rosa_aggregated_EnsemblId_subset.tsv", sep="\t")
    samples.to_csv("data/Rosa/Rosa_aggregated_samples_subset.tsv", sep="\t")
    
createSmallerAggregatedData()


# In[7]:


# Show conda environment which was used to run this notebook, as well as versions of key packages
get_ipython().system('conda info')
print("pandas version", pandas.__version__)
print("scanpy version", scanpy.__version__)

