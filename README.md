# 2018_Bacterial_Pipeline_riboseq
Developed by Fuad Mohammad
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
    
   - Bowtie_indices = index files for bowtie
 
Running the pipeline:

  - First run Github-Ribo_Density.ipynb to create 3' aligned density files. Make sure input FASTQ files are in the right directory. Specify the name of the FASTQ in the notebook or add the library to Library_names.csv. Filtering the FASTQ is done using skewer. Alignment is done using Bowtie using index files provided, or you can make your own indices. Reads mapped are then condensed to the 3' end and saved as a dictionary. 
  - Second, run Github-Ribo_Analysis.ipynb to analyze ribosome footprints at the genomic level using gene averaged information. The annotation file (Coli.GFF) is converted into a dictionary to speed up analysis, and is used to identify gene positions on the genome. 
