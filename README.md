# 2018_Bacterial_Pipeline_riboseq
  
Dependencies:
  - Jupyter notebook
  - python 2.7
  - skewer v0.2.2 
  - bowtie v0.12.7
  - BCBio
  - Biopython
  
Files to run the pipeline:
  - Github-Ribo_Density.ipynb  = pipeline to process FASTQ, align and create mapped read density files 
  - Github-Ribo_Analysis.ipynb = pipeline to analyze and visualize data
  
  - ribo_main.py     = python script containing functions for Github-Ribo_Density.ipynb
  - ribo_analysis.py = python script containing functions for Github-Ribo_Analysis.ipynb
  - ribo_util.py     = python script containing common functions for file management and data analysis
  - ribo_util.py     = python script containing plotting functions
  
  - Library_names.csv = template .csv file that can be used to keep track of libraries being analyzed. 
                        Pipeline will look for library names from here
                        
  - Annotation folder = contains annotation information for E. coli. not all of these are used.
    - bad_genes.csv      = pseudogenes, genes with homologues, insertion elements, etc. 
    - annotation.csv     = used to categorize genes by function (toxin/antitoxin, ribosomal protein, etc)
    - localization.csv   = identify genes interacting with SRP for membrane localization
    - protein_coding.csv = genes that code for proteins
    - Coli.gff           = E. coli ref genome NC_000913.2 annotation file
    - Coli_dict          = python dictionary of annotation information (gene name, seq, start, stop, etc)
  
