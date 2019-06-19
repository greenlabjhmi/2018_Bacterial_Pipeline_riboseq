'''
    Ribosome density analysis script, to analyze features of 3' or 5' mapped ribosome density
    Copyright (C) 2019  Fuad Mohammad, fuadm424@gmail.com


    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
'''

import csv
from Bio.Seq import Seq
import os 
import textwrap
import pandas as pd
from datetime import datetime
import seaborn as sns
import matplotlib.pyplot as plt
import itertools
import cPickle as pickle
import numpy as np
import ribo_util
from IPython.display import display

''' Analysis uses pickled dictionaries of densities separated by read size {length: [0,1,3,10,0...]}
    Analysis also uses pickled GFF dictionaries, Keys = [Alias, Start, Stop, Strand, Sequence]
        - Sequence has 50nt flanking the annotated CDS 
    Data for analyis is output with all settings variables as part of its name.


Table of Contents:

    each function used for analysis has a wrapper function allowing for parallel processing of many samples. 
    Make sure to consider RAM limitations under parallel processing. 

    -- SINGLE LIBRARY ANALYSIS:
    
        - readQC
        - avggenes
        - frame
        - pausescore
        - asymmetry
        - genelist
        
    -- LIBRARY COMPARISON ANALYSIS:
    
        - TE
 
        
'''



def filelen(fName):
    
    '''
    counts number of lines in a file, useful for readQC of FastQ file
    '''
    
    with open(fName) as f:
        for i, l in enumerate(f):
            pass
    
    f.close()
    return i + 1

def readQC(inputs, paths_in, paths_out):
    
    '''wrapper function for run_readQC'''
    
    
    files     = inputs['files']
    minlength = inputs['minlength']
    maxlength = inputs['maxlength']
    threads   = inputs['threads'] 
    run       = inputs['run_readQC']
    arguments = []
    
    if not files:
        print("There are no files")
        return
    
    for fname in files:
                
        path_chr = paths_out['path_chr'] + fname + '_match.fastq'
        path_analysis = paths_out['path_analysis'] + fname + '/readQC/'
        if not os.path.exists(path_analysis):
            os.makedirs(path_analysis)
            
        if not os.path.exists(path_chr):
            print "ERROR:" + fname + " has not been aligned"
            return inputs
        
        if not run == 'yes':
            if not os.path.exists(path_analysis + 'read_distribution'):
                print "ERROR:" + fname + " has not been QC analyzed"
                continue
            else: 
                print fname + " has been QC analyzed"
                continue
                
        argument = [fname, path_chr, path_analysis, minlength, maxlength]
        arguments.append(argument)
    
    print "\n\tStarted readQC at " + str(datetime.now())    
    ribo_util.multiprocess(run_readQC, arguments, threads)
    print "\tFinished readQC at " + str(datetime.now())
    
    return 


def run_readQC(fname, fastq_in, pickle_out, minlength, maxlength):
    
    '''
    input:  chromosome aligned fastq
    output: nucleotide composition per position (from the 3' end)
            read size distribution
    '''
    
    '''Get file length'''
    
    file_length = filelen(fastq_in) 
    
    '''Open file'''

    f = open(fastq_in)
    
    
    '''Define Variables and datactructures'''
    
    lengthindex = range(minlength, maxlength+1)
    
    length_data = {length : 0  for length in lengthindex}
    length_dist = {length : 0  for length in lengthindex}
    
    # add information to dictionaries: each nucleotide has a dict:
    # G_list{length : [ position = number of reads w/ G in this position]} 
    # each lenght has a list with positional composition of reads, starting at the 3' end
        
    G_data = {length : [0]*(maxlength+1) for length in lengthindex}
    C_data = {length : [0]*(maxlength+1) for length in lengthindex}
    A_data = {length : [0]*(maxlength+1) for length in lengthindex}
    U_data = {length : [0]*(maxlength+1) for length in lengthindex}

    ''' go through fastq file line by line to get info for each read: '''
    
    # each read has 4 lines of info, so divide by 4:
    
    for lines in range(0,file_length / 4):
        # each line has unique information:
        
        Identifier  = f.readline().rstrip("\n")
        Sequence    = f.readline().rstrip("\n")
        QIdentifier = f.readline().rstrip("\n")
        PHRED       = f.readline().rstrip("\n")
        
        '''count reads for each readlength'''
        
        seq_len = len(Sequence) 
        length_data[seq_len] += 1       
        
        '''calculate nucleotide composition at each position from 3' end, for each length'''
        
        # Deconstruct each read by individual nucleotides from the 3' end, noting how long the read is. 
        # Then count, adding count to dictionary with matching nucleotide identity and into the key 
        # matching the readlenght, and into the list position (index) cooresponting to the nucleotide's 
        # position from the 3' end.
        
        for position in range(1, maxlength+1):
            if seq_len >= position:
                if Sequence[-position] == 'G':
                    G_data[seq_len][position-1] += 1
                elif Sequence[-position] == 'C':
                    C_data[seq_len][position-1] += 1
                elif Sequence[-position] == 'A':
                    A_data[seq_len][position-1] += 1
                elif Sequence[-position] == 'T':
                    U_data[seq_len][position-1] += 1

    f.close()
    
    # transform absolute count to fraction of total in a given lenght of reads
    
    for length in lengthindex:  
        if length_data[length] >=1:
            T = length_data[length]    # T = total number of nucleotides of a length sequenced 
        else:
            T = 1
            
        # divide positional read count by total to get fraction of total
        
        for position in range(0, maxlength):
                
            g = float(G_data[length][position]) / float(T) *100
            c = float(C_data[length][position]) / float(T) *100 
            a = float(A_data[length][position]) / float(T) *100 
            u = float(U_data[length][position]) / float(T) *100 
            
            # replace absolute counts with fractions calculated above:
            
            G_data[length][position] = g
            C_data[length][position] = c
            A_data[length][position] = a
            U_data[length][position] = u
    
    '''get fractional read size distribution'''
    
    # calculate total reads in library
    
    sum_reads = sum(length_data.values())
    sum_reads = float(sum_reads)
    
    # calculate read length distribution: 
    
    for length in length_dist.keys():
        sum_lenght = float(length_data[length])
        read_fraction = sum_lenght / float(sum_reads) * 100
        length_dist[length] = read_fraction
        
    '''Save data in analysis folder'''
    
    ribo_util.makePickle(G_data, pickle_out + 'comp_G')
    ribo_util.makePickle(C_data, pickle_out + 'comp_C')
    ribo_util.makePickle(A_data, pickle_out + 'comp_A')
    ribo_util.makePickle(U_data, pickle_out + 'comp_U')
    ribo_util.makePickle(length_dist, pickle_out + 'read_distribution')

    f.close()
    
    return



def run_APE_site_codons(fname, fastq_in, pickle_out):
    
    '''
    go through each read and identify the codons in the A, P , and E sites. 
    
    input:  chromosome aligned fastq
    output: reads translated into aa
    '''
    
    '''Get file length'''
    
    file_length = filelen(fastq_in) 
    
    '''Open file'''

    f = open(fastq_in)
    
    
    '''Define Variables and datactructures'''
    
    APE_aa_data = []
    
    ''' go through fastq file line by line to get info for each read: '''
    
    # each read has 4 lines of info, so divide by 4:
    
    for lines in range(0,file_length / 4):
        Identifier  = f.readline().rstrip("\n")
        Sequence    = f.readline().rstrip("\n")
        QIdentifier = f.readline().rstrip("\n")
        PHRED       = f.readline().rstrip("\n")
        
        if len(Sequence) < 24:
            continue
                                     
        codons_seq = Sequence[-12: -21: -1]
        aa_code, codon_code = ribo_util.get_genetic_code()
        
        if 'N' in codons_seq:
            continue
            
        APE_codon = [codons_seq[i:i+3] for i in range(0, 9, 3)]
        
       

        APE_aa = ''
        for codon in APE_codon:
            aa = codon_code[codon]
            APE_aa += aa
        
        APE_aa_data.append(APE_aa)
    for x in range(0, 5):
        print APE_aa_data[x]

    f.close()
    
    return


def APE_site_codons(inputs, paths_in, paths_out):
    
    '''wrapper function for run_APE_site_codons'''
    
    
    files     = inputs['files']
    
    threads   = inputs['threads'] 
    run       = 'yes'
    arguments = []
    
    if not files:
        print("There are no files")
        return
    
    for fname in files:
                
        path_chr = paths_out['path_chr'] + fname + '_match.fastq'
        path_analysis = paths_out['path_analysis'] + fname + '/readQC/'
        if not os.path.exists(path_analysis):
            os.makedirs(path_analysis)
            
        if not os.path.exists(path_chr):
            print "ERROR:" + fname + " has not been aligned"
            return inputs
        
        if not run == 'yes':
            if not os.path.exists(path_analysis + 'read_distribution'):
                print "ERROR:" + fname + " has not been QC analyzed"
                continue
            else: 
                print fname + " has been QC analyzed"
                continue
                
        argument = [fname, path_chr, path_analysis]
        arguments.append(argument)
    
    print "\n\tStarted readQC at " + str(datetime.now())    
    ribo_util.multiprocess(run_APE_site_codons, arguments, threads)
    print "\tFinished readQC at " + str(datetime.now())
    
    return 


def avggenes(inputs, paths_out, settings, gff_dict, plus_dict, minus_dict): 
    
    '''wrapper function for run_avggenes'''
    
    files     = inputs['files']
    threads   = inputs['threads']
    multi     = inputs['multiprocess']
    arguments = []
    
    # check if file list has items
    
    if not files:
        print("There are no files")
        return
    
    print "Started avggenes at " + str(datetime.now())
    
    # feed each file into function, paralleliing when desired

    for fname in files:
        
        #define output paths:
        
        path_analysis = paths_out['path_analysis'] + fname + '/avggenes/'
        if not os.path.exists(path_analysis):
            os.makedirs(path_analysis)
        path_start = path_analysis + 'avg_start_'
        path_stop  = path_analysis + 'avg_stop_'
        plus       = plus_dict[fname]
        minus      = minus_dict[fname] 
        
        # run function: if in series then begin running. Else, make arguemnt list for each fname
        
        if not multi == 'yes':
            run_avggenes(fname, settings, plus, minus, gff_dict, path_start, path_stop)
        else:    
            argument = [fname, settings, plus, minus, gff_dict, path_start, path_stop]
            arguments.append(argument)
    
    # feed argument list into multiprocessing if multi == yes
    
    if multi == 'yes':
        ribo_util.multiprocess(run_avggenes, arguments, threads)
        
    print "Finished avggenes at " + str(datetime.now())

    return

def run_avggenes(fname, settings, plus, minus, gff, path_start, path_stop):
    
    '''
    input:  Density dictionary, GFF dictionary
    Output: Average occupancy over all genes, aligned to either start or stop positions: 
    can be represented as a summation of reads, '''
    
    
    '''settings'''
    
    # define setting variables:
    minlength      = settings['minlength']
    maxlength      = settings['maxlength']
    length_in_ORF  = settings['length_in_ORF']
    length_out_ORF = settings['length_out_ORF']
    density_type   = settings['density_type'] 
    alignment      = settings['alignment']
    next_gene      = settings['next_gene']
    equal_weight   = settings['equal_weight']
    threshold      = settings['threshold']
    
    window_length  = length_in_ORF + length_out_ORF
    positionindex  = range(0, window_length) 
    lengthindex    = range(minlength, maxlength+1)
    
    # filename annotation to maintain separate figures
    minlength_1      = str(minlength)     +'_'
    maxlength_1      = str(maxlength)     +'_'
    length_in_ORF_1  = str(length_in_ORF) +'_'
    length_out_ORF_1 = str(length_out_ORF)+'_'
    density_type_1   = density_type       +'_'
    next_gene_1      = str(next_gene)     +'_'
    equal_weight_1   = equal_weight       +'_'
    threshold_1      = str(threshold)     +'_'
    
    name_settings = length_in_ORF_1+length_out_ORF_1+next_gene_1+threshold_1
    name_settings += density_type_1+equal_weight_1+minlength_1+maxlength_1
    
    # define density variables:
    density_plus  = plus
    density_minus = minus
    
    if density_type == 'rpm':
        #convert read to rpm
        density_plus, density_minus = ribo_util.get_density_rpm(plus, minus, minlength, maxlength)
        
    totalcounts = ribo_util.get_allcounts(density_plus, density_minus, minlength, maxlength)
        
    # define annotation variables:    
    gff_dict = gff
    alias_list  = gff_dict['Alias'] 
    strand_list = gff_dict['Strand'] 
    start_list  = gff_dict['Start'] 
    stop_list   = gff_dict['Stop'] 
    #type_list   = gff_dict['Type']
    
    
    '''datastructures:'''
    
    # for heatmap - each dict has read counts at each position separated by length keys    
    averagegene_start = {length : [0]*(window_length) for length in lengthindex}
    averagegene_stop  = {length : [0]*(window_length) for length in lengthindex}
    
    # for avggene summary - dictionary with position as keys 
    start_all = {position : 0 for position in positionindex}
    stop_all  = {position : 0 for position in positionindex}
    
    # save CDS vs 3'UTR density ratio for every gene that passes thresholds:
    utr_vs_cds = {} # {alias: ratio}
    
    # count genes excluded from data  = [count, [names of genes]]   
    excluded_genes = {}
    excluded_genes['short']     = [0, []]
    excluded_genes['threshold'] = [0, []]
    excluded_genes['too_close'] = [0, []]
    excluded_genes['type']      = [0, []]
    # count included genes
    included_genes = [0, []]
    
    
    '''Calculate gene averages and add to data structure'''
    
    for alias, start, stop, strand in itertools.izip(alias_list, start_list,stop_list, strand_list):  
        genelength = abs(start - stop)
        
        nextgene = ribo_util.nextgene(alias, gff_dict) # outputs {'distance': num, 'alias': alias}
        #prevgene = ribo_util.prevgene(alias, gff_dict) # outputs {'distance': num, 'alias': alias}
        
        # define plot start and stop window for genes in + or - strand:
        
        if strand == '+':
            start_window = start - length_out_ORF 
            stop_window  = stop - length_in_ORF 
            stop_window_max = stop + length_out_ORF 
            utr_start = stop + 20
            utr_end   = stop_window_max
            
            
        if strand == '-':
            start_window = start + length_out_ORF 
            stop_window  = stop + length_in_ORF 
            stop_window_max = stop - length_out_ORF 
            utr_start = stop - 20
            utr_end   = stop_window_max - 20
               
        # exclude genes that are too close 
            
        if alias in excluded_genes['too_close'][1]:
            excluded_genes['too_close'][0] += 1
            continue 
        
        if next_gene != 'no':
            if nextgene['distance'] < next_gene: 
                #print nextgene['distance'], alias, nextgene['alias']
                excluded_genes['too_close'][0] += 1
                excluded_genes['too_close'][1].append(alias)
                excluded_genes['too_close'][1].append(nextgene['alias'])
                continue
                
        #exclude genes that are not CDS
        '''feat_type = 'CDS'
        if feat_type != 'CDS':
            excluded_genes['type'][1].append(alias)
            excluded_genes['type'][0] += 1
            continue'''
        
        #exclude genes that are too small    
        if genelength < length_in_ORF :  
            excluded_genes['short'][0] += 1
            excluded_genes['short'][1].append(alias)
            continue 
        
        '''Get read counts'''
        
        # reads on CDS:
        gene_reads, length_reads = ribo_util.get_genecounts(start, stop, strand, density_plus, 
                                                             density_minus, minlength, maxlength)
        gene_reads = [float(i) for i in gene_reads] 
        
        # reads over window to plot:
        window_start = start_window 
        window_stop  = stop_window_max
        windowlength = abs(start_window - stop_window_max)
        
        window_reads, window_length_reads = ribo_util.get_genecounts(window_start, window_stop, strand, density_plus, 
                                                                        density_minus, minlength, maxlength)
        window_reads = [float(i) for i in window_reads] 
        
        # reads over UTR, length determined by plot settings
        utrlength = abs(utr_end - utr_start)
        
        utr_reads, utr_length_reads = ribo_util.get_genecounts(utr_start, utr_end, strand, density_plus, 
                                                                        density_minus, minlength, maxlength)
        utr_reads = [float(i) for i in utr_reads]             

        #exclude genes with low rpkm   
        if density_type == 'rpm':
            gene_rpm  = sum(gene_reads)
            gene_rpkm = gene_rpm / float(genelength) * 1000
            gene_rpc  = gene_rpm / (float(genelength)/3) 
        else: 
            gene_rpm  = float(sum(gene_reads)) / float(totalcounts) * 1000000
            gene_rpkm = gene_rpm / float(genelength) * 1000
            gene_rpc  = gene_rpm / (float(genelength)/3) 
        
        if gene_rpc < threshold or gene_rpc == 0:
                excluded_genes['threshold'][0] += 1
                excluded_genes['threshold'][1].append(alias)
                continue
                  
        # calculate UTR/CDS read ratio
        utr_norm = sum(utr_reads) / float(utrlength)
        cds_norm = sum(gene_reads) / float(genelength)
        
        if not cds_norm == 0:
            utr_cds_ratio = utr_norm / cds_norm
        else: 
            utr_cds_ratio = 0
         
        # equal weight normalization: normalize all genes to be weighted equally in average
        if equal_weight == 'yes':
            normalization = sum(window_reads) / windowlength
        else:
            normalization = 1
       
        genomelength = len(density_plus[density_plus.keys()[0]])
        
        # add gene density to datastructure:
        if strand == '+':
            
            # remove genes that are at the beginning or end of genome if problematic
            if not 0 <= start_window < genomelength:
                excluded_genes['too_close'][0] += 1
                excluded_genes['too_close'][1].append(alias)
                continue
            if not 0 <= stop_window_max < genomelength:
                excluded_genes['too_close'][0] += 1
                excluded_genes['too_close'][1].append(alias)
                continue
            
            density_dict = density_plus
            for length in lengthindex:
                for position in positionindex:
                    start_density = density_dict[length][start_window + position] 
                    start_density = float(start_density) / normalization
                    averagegene_start[length][position] += start_density
                    start_all[position] += start_density
                    
                    stop_density = density_dict[length][stop_window + position] 
                    stop_density = float(stop_density) / normalization
                    averagegene_stop[length][position] += stop_density
                    stop_all[position]  += stop_density

        if strand == '-':
            
            if not 0 <= start_window < genomelength:
                excluded_genes['too_close'][0] += 1
                excluded_genes['too_close'][1].append(alias)
                continue
            if not 0 <= stop_window_max < genomelength:
                excluded_genes['too_close'][0] += 1
                excluded_genes['too_close'][1].append(alias)
                continue
                
            density_dict = density_minus
            for length in lengthindex:
                for position in positionindex:
                    start_density = density_dict[length][start_window - position]
                    start_density = float(start_density) / normalization
                    averagegene_start[length][position] += start_density
                    start_all[position] += start_density
                    
                    stop_density = density_dict[length][stop_window - position] 
                    stop_density = float(stop_density) / normalization
                    averagegene_stop[length][position] += stop_density 
                    stop_all[position]  += stop_density 
        
        # count included genes 
        included_genes[0] += 1
        included_genes[1].append(alias)
        
        # add utr/cds to datframe:
        utr_vs_cds[alias] = utr_cds_ratio
    
    if equal_weight == 'yes':
        gene_count = float(included_genes[0])
        norm_factor = 1 / gene_count
        
        for position in positionindex:
            for length in lengthindex:
                averagegene_start[length][position] = averagegene_start[length][position] * norm_factor 
                averagegene_stop[length][position]  = averagegene_stop[length][position] * norm_factor
            
            start_all[position] = start_all[position] * norm_factor
            stop_all[position]  = stop_all[position] * norm_factor
        
    # Print report of included and excluded genes:
    excluded   = excluded_genes['too_close'][0] + excluded_genes['short'][0] + excluded_genes['threshold'][0] + excluded_genes['type'][0]
    
    totalgenes = excluded + included_genes[0]
    
    print '\tFor ' + fname + ': ' + str(included_genes[0]) + ' out of ' + str(totalgenes) + ' genes in average.'
    print '\t    - genes failing distance cutoff = ' + str(excluded_genes['too_close'][0])
    print '\t    - genes failing length cutoff   = ' + str(excluded_genes['short'][0])
    print '\t    - genes failing RPM cutoff      = ' + str(excluded_genes['threshold'][0])    
    print '\t    - non-coding genes              = ' + str(excluded_genes['type'][0]) 
    
    
    # Save pickled versions of the datastructures:
    # name_settings: inorf_outorf_nextgene_threshold_reads/rpm_weightedyes/no_minlength_maxlength
    
    ribo_util.makePickle(start_all,         path_start + name_settings + '_all', protocol=pickle.HIGHEST_PROTOCOL)              
    ribo_util.makePickle(averagegene_start, path_start + name_settings + '_HM' , protocol=pickle.HIGHEST_PROTOCOL)
    ribo_util.makePickle(stop_all,          path_stop + name_settings + '_all' , protocol=pickle.HIGHEST_PROTOCOL)
    ribo_util.makePickle(averagegene_stop,  path_stop + name_settings + '_HM'  , protocol=pickle.HIGHEST_PROTOCOL)
    ribo_util.makePickle(utr_vs_cds,        path_start + name_settings + '_UTR_CDS', protocol=pickle.HIGHEST_PROTOCOL)
    
    return

                



def run_frame(fname, settings, plus, minus, gff, path_analysis):
    
    # define setting variables:
    minlength = settings['minlength']
    maxlength = settings['maxlength']
    threshold = settings['threshold']
    
    lengthindex   = range(minlength, maxlength+1)

    # define density variables:
    density_plus  = plus
    density_minus = minus
    totalcounts   = ribo_util.get_allcounts(density_plus, density_minus, minlength, maxlength)
        
    # define annotation variables:    
    gff_dict = gff
    alias_list  = gff_dict['Alias'] 
    strand_list = gff_dict['Strand'] 
    start_list  = gff_dict['Start'] 
    stop_list   = gff_dict['Stop'] 
    
    # datastructure:
    length_frame   = {length : [0, 0, 0] for length in lengthindex}
    alias_frame    = {}
    genome_frame   = [0, 0, 0]
    genome_total   = 0
    
    # count genes excluded from data  = [count, [names of genes]]   
    excluded_genes = {}
    excluded_genes['short']     = [0, []]
    excluded_genes['threshold'] = [0, []]
    # count included genes
    included_genes = [0, []]
    
    for alias, start, stop, strand in itertools.izip(alias_list, start_list,stop_list, strand_list):  
        
        genelength = abs(start - stop) + 1   
       
        gene_reads, length_reads = ribo_util.get_genecounts(start, stop, strand, density_plus, 
                                                            density_minus, minlength, maxlength)
        gene_reads = [float(i) for i in gene_reads]         

        #exclude genes with low rpkm   
        gene_rpm  = float(sum(gene_reads)) / float(totalcounts) * 1000000
        gene_rpkm = gene_rpm / float(genelength) * 1000
        gene_rpc  = gene_rpm / float(genelength) * 3
        
        if gene_rpkm < threshold or gene_rpkm == 0:
            excluded_genes['threshold'][0] += 1
            excluded_genes['threshold'][1].append(alias)
            continue
        
        if genelength < 50:
            excluded_genes['short'][0] += 1
            excluded_genes['short'][1].append(alias)
            continue
            
        if strand == '+':
            density_dict = density_plus
            startframe   = start + 27 
            stopframe    = stop  - 12
            
            start_index_0 = startframe
            start_index_1 = startframe + 1
            start_index_2 = startframe + 2
            
            stop_index = stopframe + 1
            period = 3
            
        elif strand == '-':
            density_dict = density_minus
            startframe   = start - 27
            stopframe    = stop  + 12 
            
            start_index_0 = startframe
            start_index_1 = startframe - 1
            start_index_2 = startframe - 2
            
            stop_index = stopframe - 1
            period = -3
        
        alias_frame[alias] = [0, 0, 0]
        
        for length in lengthindex:
            frame_0 = density_dict[length][start_index_0: stop_index + 1: period]
            frame_1 = density_dict[length][start_index_1: stop_index + 1: period]
            frame_2 = density_dict[length][start_index_2: stop_index + 1: period]

            frame_0 = float(sum(frame_0))
            frame_1 = float(sum(frame_1))
            frame_2 = float(sum(frame_2))
            
            length_frame[length][0] += frame_0
            length_frame[length][1] += frame_1
            length_frame[length][2] += frame_2

            alias_frame[alias][0] += frame_0
            alias_frame[alias][1] += frame_1
            alias_frame[alias][2] += frame_2

            genome_frame[0] += frame_0
            genome_frame[1] += frame_1
            genome_frame[2] += frame_2
                
        alias_total   = sum(alias_frame[alias])
        genome_total += alias_total
            
        if alias_total <= 0:
            alias_total = 1
                
        alias_frame[alias][0] = 100 * alias_frame[alias][0] / alias_total
        alias_frame[alias][1] = 100 * alias_frame[alias][1] / alias_total
        alias_frame[alias][2] = 100 * alias_frame[alias][2] / alias_total
               
    for length in lengthindex:
        length_total = sum(length_frame[length])
        if length_total <= 0:
            length_total = 1
                        
        length_frame[length][0] = 100 * length_frame[length][0] / length_total
        length_frame[length][1] = 100 * length_frame[length][1] / length_total
        length_frame[length][2] = 100 * length_frame[length][2] / length_total
    
    genome_frame[0] = 100 * genome_frame[0] / genome_total
    genome_frame[1] = 100 * genome_frame[1] / genome_total
    genome_frame[2] = 100 * genome_frame[2] / genome_total

    # filename annotation to maintain separate figures
    minlength_1      = str(minlength)     +'_'
    maxlength_1      = str(maxlength)     +'_'
    threshold_1      = str(threshold)     +'_'
    
    name_settings = minlength_1+maxlength_1+threshold_1
    
    ribo_util.makePickle(alias_frame,  path_analysis + name_settings + 'frame_alias' , protocol=pickle.HIGHEST_PROTOCOL) 
    ribo_util.makePickle(length_frame, path_analysis + name_settings + 'frame_length', protocol=pickle.HIGHEST_PROTOCOL)    
    ribo_util.makePickle(genome_frame, path_analysis + name_settings + 'frame_genome', protocol=pickle.HIGHEST_PROTOCOL)
    
    return

        
def frame(inputs, paths_out, settings, gff_dict, plus_dict, minus_dict): 
    files     = inputs['files']
    threads   = inputs['threads']
    multi     = inputs['multiprocess']
    arguments = []
    
    if not files:
        print("There are no files")
        return
    
    print "Started frame analysis at " + str(datetime.now())

    for fname in files:
        path_analysis = paths_out['path_analysis'] + fname + '/frame/'
        if not os.path.exists(path_analysis):
            os.makedirs(path_analysis)
        plus  = plus_dict[fname]
        minus = minus_dict[fname] 
        
        if not multi == 'yes':
            run_frame(fname, settings, plus, minus, gff_dict, path_analysis)
        else:     
            argument = [fname, settings, plus, minus, gff_dict, path_analysis]
            arguments.append(argument)
    
    if multi == 'yes':
        ribo_util.multiprocess(run_frame, arguments, threads)
    print "Finished frame analysis at " + str(datetime.now())
    
    return
    

def run_pausescore(fname, settings, plus, minus, gff, path_pausescore):
    
    '''define variables'''
    
    minlength   = settings['minlength']
    maxlength   = settings['maxlength']
    alignment   = settings['alignment']
    threshold   = settings['threshold']


    lengthindex = range(minlength, maxlength + 1)
    
    A_site = settings['A_site shift']
    P_site = A_site - 3
    E_site = A_site - 6
    two_site = A_site + 6
    one_site = A_site + 3
    mone_site = A_site - 9
    mtwo_site = A_site - 12
    
    
    frameshift      = settings['frameshift']

    plot_upstream   = settings['plot_upstream'] / 3 * 3        #change window to interval of 3
    plot_downstream = settings['plot_downstream'] / 3 * 3
    
    next_codon = settings['next_codon']
    
    
    # define plot length
    plotlength = plot_upstream + plot_downstream + 1
    positionindex = range(0, plotlength)
    
    # load density files
    density_plus  = plus
    density_minus = minus 
    
    # load annotation
    gff_dict = gff
    
    alias_list  = gff_dict['Alias'] 
    strand_list = gff_dict['Strand'] 
    start_list  = gff_dict['Start'] 
    stop_list   = gff_dict['Stop'] 
    seq_list    = gff_dict['Sequence'] 
 
    gff_extra  = settings['gff_extra']   # extra nucleotides in gff_dict sequence (UTR sequence)
    start_trim = settings['start_trim'] / 3 * 3
    stop_trim  = settings['stop_trim'] / 3 * 3   
    
    
    minlength_1       = str(minlength)      +'_'
    maxlength_1       = str(maxlength)      +'_'
    plot_upstream_1   = str(plot_upstream)  +'_'
    plot_downstream_1 = str(plot_downstream)+'_'
    start_trim_1      = str(start_trim)     +'_'
    stop_trim_1       = str(stop_trim)      +'_'
    frameshift_1      = str(frameshift)     +'_'
    a_site_1          = str(A_site)         +'_'
    
    name_settings = minlength_1+maxlength_1+plot_upstream_1+plot_downstream_1
    name_settings += start_trim_1+stop_trim_1+frameshift_1
        
    # import genetic code
    aa_code, codon_code = ribo_util.get_genetic_code()
    
    
    '''output data structure'''
    
    # aa/codon avgplot      = { aa: {length: [values]}} 
    # aa/codon count        = { aa: N}
    # aa/codon score        = { aa: [aa_list], A site score: [A_list]...
                                
    aa_avgplot  = {}
    aa_count    = {}
    aa_score    = {}
    aa_score_df = {}
    aa_list     = []
    aa_A_list   = []
    aa_P_list   = []
    aa_E_list   = []
    aa_1_list   = []
    aa_2_list   = []
    aa_m1_list   = []
    aa_m2_list   = []
    
    codon_avgplot  = {}
    codon_count    = {}
    codon_score    = {}
    codon_score_df = {}
    codon_list     = []
    codon_A_list   = []
    codon_P_list   = []
    codon_E_list   = []
    codon_1_list   = []
    codon_2_list   = []
    codon_m1_list   = []
    codon_m2_list   = []
    
    # create empty datastructures to store density info:
    for aa in aa_code.keys():
        aa_avgplot[aa] = {length : [0]*(plotlength) for length in lengthindex}
        aa_count[aa]   = 0
        aa_score[aa]   = [0, 0, 0, 0, 0, 0, 0] # [Asite, P site, E site] 
        
    for codon in codon_code.keys():
        codon_avgplot[codon] = {length : [0]*(plotlength) for length in lengthindex}
        codon_count[codon]   = 0
        codon_score[codon]   = [0, 0, 0, 0, 0, 0, 0] # [Asite, P site, E site] 
    
    '''genes in data''' 
    
    # count genes excluded from data  = [count, [names of genes]]   
    excluded_genes = {}
    excluded_genes['short']         = [0, []]
    excluded_genes['low_density']   = [0, []]
    excluded_genes['not_divisible'] = [0, []]
    # count included genes
    included_genes = [0, []]
        
    '''iterate through every annotated gene to get codon density info:'''
    
    for alias, start, stop, strand, sequence in itertools.izip(alias_list, start_list,stop_list, strand_list, seq_list):

        
        ''' define start and stop positions for codons to analyze:
        
        # = codon 
        #################################### = GFF sequence  (50 extra nt)
           ##############################    = AUG to UGA
             ##########################      = density to analyze : remove start and stop peaks
                 ##################          = codons to analyze : remove plot window
        
        '''
        
        # First, define density without start and stop peaks:
        if strand == '+':
            density_dict  = density_plus
            density_start = start + start_trim + frameshift
            density_stop  = stop  - stop_trim + frameshift
            
            period = 1
            
        elif strand == '-':
            density_dict  = density_minus
            density_start = start - start_trim - frameshift
            density_stop  = stop  + stop_trim - 3 + frameshift
            
            period = -1
        
        # GFF seq has 50 extra nucleotides, so remove:
        # Also remove several codons from start and stop positions:
        
        codon_seq_start = gff_extra + start_trim + plot_upstream + frameshift
        codon_seq_stop  = -gff_extra - stop_trim - plot_downstream + frameshift
            
        seq       = sequence[codon_seq_start : codon_seq_stop]
        seqlength = len(seq)
        
        # exclude genes that are not divisable by 3
        if seqlength % 3 != 0:
            excluded_genes['not_divisible'][0] += 1
            excluded_genes['not_divisible'][1].append(alias)
            continue
        
        # exclude genes shorter than plot
        if seqlength < plotlength + 1:
            excluded_genes['short'][0] += 1
            excluded_genes['short'][1].append(alias)
            continue
            
        elif codon_seq_start - codon_seq_stop > len(sequence):
            excluded_genes['short'][0] += 1
            excluded_genes['short'][1].append(alias)
            continue
    
        # make a list of codons in the sequence:
        codons_seq = [seq[i:i+3] for i in range(0, seqlength, 3)]
        aa_seq     = [codon_code[codon] for codon in codons_seq]
        
        #make empty density dict for the gene
        genelength    = abs(density_start - density_stop) + 1
        gene_density  = {}
        total_density = [0] * genelength
        
        # fill density dict with density info
        # gives density encompassed by seq, plus extra defined by plotlength
        for length in lengthindex:
            length_density       = density_dict[length][density_start: density_stop: period]
            length_density_float = [float(i) for i in length_density]
            gene_density[length] = length_density
            total_density        = [x + y for x, y in itertools.izip(total_density, length_density)]
            
        gene_avgreads = float(sum(total_density)) / float(genelength)

            
        # normalize density by total gene density
        if gene_avgreads * 3 < threshold:
            excluded_genes['low_density'][0] += 1
            excluded_genes['low_density'][1].append(alias)
            continue
            
        else: 
            relative_density = {}
            
            for length in lengthindex:
                relative_density[length] = [reads / gene_avgreads for reads in gene_density[length]]
        
        # add data to dataframe
        codon_index = 0
        for codon in codons_seq:
            codon_position = (codon_index * 3)
            aa = aa_seq[codon_index]
            codon_index        += 1
            
            '''if codon_index > 50:
                continue'''
            
            codon_count[codon] += 1
            aa_count[aa]       += 1
            
            for length in lengthindex:
                for position in positionindex:
                    
                    density_position = codon_position + position
                    density          = relative_density[length][density_position]
                    
                    codon_avgplot[codon][length][position] += density
                    aa_avgplot[aa][length][position]       += density
                    
        included_genes[0] += 1
        included_genes[1].append(alias)
    
    #divide data by total instances of the codon or aa
    for codon in codon_avgplot.keys():
        for length in lengthindex:
            for position in positionindex:
                if codon_count[codon] == 0:
                    continue 
                else: 
                    codon_avgplot[codon][length][position] = codon_avgplot[codon][length][position] / codon_count[codon]
    
    for aa in aa_avgplot.keys():
        for length in lengthindex:
            for position in positionindex:
                if aa_count[aa] == 0:
                    continue
                else: 
                    aa_avgplot[aa][length][position] = aa_avgplot[aa][length][position] / aa_count[aa]
                
    '''Convert data for plotting and csv'''
    
    codon_data     = {}
    codon_data_sum = {}
    aa_data     = {}
    aa_data_sum = {}
    
    A_site_shift = plot_upstream - A_site 
    P_site_shift = plot_upstream - P_site 
    E_site_shift = plot_upstream - E_site
    one_site_shift = plot_upstream - one_site 
    two_site_shift = plot_upstream - two_site 
    mone_site_shift = plot_upstream - mone_site 
    mtwo_site_shift = plot_upstream - mtwo_site 
    
    for codon in codon_avgplot.keys():
        df = pd.DataFrame(codon_avgplot[codon])
        df = df.T
        df = df.reindex(index=df.index[::-1])

        plot_range  = len(df.columns)
        shift_start = - plot_upstream
        shift_stop  = plot_range - plot_upstream

        df.columns = pd.RangeIndex(start = shift_start, stop = shift_stop, step = 1)

        codon_data[codon]     = df
        codon_data_sum[codon] = df.sum(0)
        
        codon_score[codon][0] = sum(codon_data_sum[codon][A_site_shift: A_site_shift + 3:1])/3
        codon_score[codon][1] = sum(codon_data_sum[codon][P_site_shift: P_site_shift + 3:1])/3                            
        codon_score[codon][2] = sum(codon_data_sum[codon][E_site_shift: E_site_shift + 3:1])/3
        codon_score[codon][3] = sum(codon_data_sum[codon][one_site_shift: one_site_shift + 3:1])/3
        codon_score[codon][4] = sum(codon_data_sum[codon][two_site_shift: two_site_shift + 3:1])/3                            
        codon_score[codon][5] = sum(codon_data_sum[codon][mone_site_shift: mone_site_shift + 3:1])/3
        codon_score[codon][6] = sum(codon_data_sum[codon][mtwo_site_shift: mtwo_site_shift + 3:1])/3
        
        codon_list.append(codon)
        codon_A_list.append(codon_score[codon][0])
        codon_P_list.append(codon_score[codon][1])
        codon_E_list.append(codon_score[codon][2])
        codon_1_list.append(codon_score[codon][3])
        codon_2_list.append(codon_score[codon][4])
        codon_m1_list.append(codon_score[codon][5])
        codon_m2_list.append(codon_score[codon][6])

    for aa in aa_avgplot.keys():
        df = pd.DataFrame(aa_avgplot[aa])
        df = df.T
        df = df.reindex(index=df.index[::-1])

        plot_range  = len(df.columns)
        shift_start = - plot_upstream
        shift_stop  = plot_range - plot_upstream

        df.columns = pd.RangeIndex(start = shift_start, stop = shift_stop, step = 1)

        aa_data[aa]     = df
        aa_data_sum[aa] = df.sum(0) 
                                    
        aa_score[aa][0] = sum(aa_data_sum[aa][A_site_shift: A_site_shift + 3:1])/3    
        aa_score[aa][1] = sum(aa_data_sum[aa][P_site_shift: P_site_shift + 3:1])/3
        aa_score[aa][2] = sum(aa_data_sum[aa][E_site_shift: E_site_shift + 3:1])/3
        aa_score[aa][3] = sum(aa_data_sum[aa][one_site_shift: one_site_shift + 3:1])/3    
        aa_score[aa][4] = sum(aa_data_sum[aa][two_site_shift: two_site_shift + 3:1])/3
        aa_score[aa][5] = sum(aa_data_sum[aa][mone_site_shift: mone_site_shift + 3:1])/3 
        aa_score[aa][6] = sum(aa_data_sum[aa][mtwo_site_shift: mtwo_site_shift + 3:1])/3 

        aa_list.append(aa)
        aa_A_list.append(aa_score[aa][0])
        aa_P_list.append(aa_score[aa][1])
        aa_E_list.append(aa_score[aa][2])
        aa_1_list.append(aa_score[aa][3])
        aa_2_list.append(aa_score[aa][4])
        aa_m1_list.append(aa_score[aa][5])
        aa_m2_list.append(aa_score[aa][6])
        
    codon_score_df['Codon']  = codon_list
    codon_score_df['A_site'] = codon_A_list
    codon_score_df['P_site'] = codon_P_list
    codon_score_df['E_site'] = codon_E_list
    codon_score_df['1_site'] = codon_1_list
    codon_score_df['2_site'] = codon_2_list
    codon_score_df['-1_site'] = codon_m1_list
    codon_score_df['-2_site'] = codon_m2_list

    aa_score_df['Amino Acid'] = aa_list
    aa_score_df['A_site']     = aa_A_list
    aa_score_df['P_site']     = aa_P_list
    aa_score_df['E_site']     = aa_E_list
    aa_score_df['1_site']     = aa_1_list
    aa_score_df['2_site']     = aa_2_list
    aa_score_df['-1_site']     = aa_m1_list
    aa_score_df['-2_site']     = aa_m2_list

    codon_df = pd.DataFrame(codon_score_df)
    codon_df.to_csv(path_pausescore + 'codon_scores.csv')
    
    aa_df = pd.DataFrame(aa_score_df)
    aa_df.to_csv(path_pausescore + 'aa_scores.csv')
    
    ribo_util.makePickle(codon_data,     path_pausescore +'codon_HM_data'  + name_settings , protocol=pickle.HIGHEST_PROTOCOL) 
    ribo_util.makePickle(codon_data_sum, path_pausescore +'codon_plot_data'+ name_settings , protocol=pickle.HIGHEST_PROTOCOL)  
    ribo_util.makePickle(codon_score_df, path_pausescore +'codon_scores'   + name_settings , protocol=pickle.HIGHEST_PROTOCOL)
    ribo_util.makePickle(aa_data,        path_pausescore +'aa_HM_data'     + name_settings , protocol=pickle.HIGHEST_PROTOCOL) 
    ribo_util.makePickle(aa_data_sum,    path_pausescore +'aa_plot_data'   + name_settings , protocol=pickle.HIGHEST_PROTOCOL)
    ribo_util.makePickle(aa_score_df,    path_pausescore +'aa_scores'      + name_settings , protocol=pickle.HIGHEST_PROTOCOL)
    
    #print excluded_genes['short'][0]    
    #print excluded_genes['low_density'][0]
    #print excluded_genes['not_divisible'][0]
    #print included_genes[0]
    return 


def pausescore(inputs, paths_out, settings, gff_dict, plus_dict, minus_dict):
    
    files     = inputs['files']
    threads   = inputs['threads'] 
    multi     = inputs['multiprocess']
    arguments = []
    
    if not files:
        print("There are no files")
        return
    
    print "Started pause score analysis at " + str(datetime.now())

    for fname in files:
        path_pausescore = paths_out['path_analysis'] + fname + '/pause_score/'
        if not os.path.exists(path_pausescore):
            os.makedirs(path_pausescore)
        plus  = plus_dict[fname]
        minus = minus_dict[fname] 
       
        if not multi == 'yes':
            run_pausescore(fname, settings, plus, minus, gff_dict, path_pausescore)
        else:     
            argument = [fname, settings, plus, minus, gff_dict, path_pausescore]
            arguments.append(argument)
    
    if multi == 'yes':
        ribo_util.multiprocess(run_pausescore, arguments, threads)
    
    print "Finished pause score analysis at " + str(datetime.now())
    
    return



                
def run_asymmetry(fname, settings, plus, minus, gff, path_analysis): 
    
    minlength    = settings['minlength']
    maxlength    = settings['maxlength']
    threshold    = settings['threshold']
    subgroup_csv = settings['subgroup']
    
    minlength_1       = str(minlength)      +'_'
    maxlength_1       = str(maxlength)      +'_'
    
    name_settings = minlength_1+maxlength_1
    
    density_plus  = plus
    density_minus = minus 
    gff_dict      = gff
    
    if not subgroup_csv == 'none':
        subgroup       = pd.read_csv(subgroup_csv)
        subgroup       = subgroup.to_dict(orient='list')
        subgroup_alias = subgroup['Alias']
    else:
        subgroup_alias = ['none']

    lengthindex = range(minlength, maxlength+1)
    
    alias_list  = gff_dict['Alias'] 
    strand_list = gff_dict['Strand'] 
    start_list  = gff_dict['Start'] 
    stop_list   = gff_dict['Stop'] 
    
    asymmetry_alias      = []
    asymmetry_subgroup   = []
    asymmetry_score      = []
    asymmetry_reads      = []
    asymmetry_genelength = []
    asymmetry_dict  = {}
    genes_excluded  = []
       
    for alias, start, stop, strand in itertools.izip(alias_list, start_list,stop_list, strand_list):  
        
        genelength = abs(start - stop) + 1
             
        if genelength < 80: 
            genes_excluded.append(alias)
            continue
            
        mid_gene = genelength / 2
        
        if strand == '+':
            density_dict = density_plus
            midpoint = start + mid_gene 
            startp = start + 27
            stopp  = stop  - 12
            period = 1
            
        elif strand == '-':
            density_dict = density_minus
            midpoint = start - mid_gene 
            startp = start - 27
            stopp  = stop  + 12
            period = -1
        
        first_half  = 0 
        second_half = 0 
        
        for length in lengthindex:
            
            first  = density_dict[length][startp: midpoint + 1: period]
            second = density_dict[length][midpoint: stopp  + 1: period]
            genelength = len(density_dict[length][start: stop + 1: period])
            
            first  = float(sum(first))
            second = float(sum(second))
            
            first_half  += first
            second_half += second
            
        if first_half == 0:
            genes_excluded.append(alias)
            continue
            
        reads =  second_half + first_half    
        
        asymmetry = second_half / first_half
        asymmetry = np.log2(asymmetry)
        
        if reads < 150:
            genes_excluded.append(alias)
            continue
                                
        if asymmetry < -3 or asymmetry > 3:
            genes_excluded.append(alias)
            continue    
        
        reads = np.log10(reads)
        genelength = np.log10(genelength)
        
        if alias in subgroup_alias:
            asymmetry_subgroup.append('Subgroup')                          
        else:
            asymmetry_subgroup.append('Other')
        
        asymmetry_alias.append(alias)
        asymmetry_score.append(asymmetry)
        asymmetry_reads.append(reads)
        asymmetry_genelength.append(genelength)
    
    asymmetry_dict['Alias']      = asymmetry_alias
    asymmetry_dict['Score']      = asymmetry_score
    asymmetry_dict['Reads']      = asymmetry_reads
    asymmetry_dict['GeneLength'] = asymmetry_genelength
    asymmetry_dict['Subgroup']   = asymmetry_subgroup
 
    ribo_util.makePickle(asymmetry_dict, path_analysis + 'asymmetry', protocol=pickle.HIGHEST_PROTOCOL) 
                
            
def asymmetry(inputs, paths_out, settings, gff_dict, plus_dict, minus_dict):
    
    files     = inputs['files']
    threads   = inputs['threads'] 
    arguments = []
    
    if not files:
        print("There are no files")
        return
    
    print "Started asymmetry analysis at " + str(datetime.now())

    for fname in files:
        path_analysis = paths_out['path_analysis'] + fname + '/'
        if not os.path.exists(path_analysis):
            os.makedirs(path_analysis)
        plus  = plus_dict[fname]
        minus = minus_dict[fname] 
        
        run_asymmetry(fname, settings, plus, minus, gff_dict, path_analysis)
        argument = [fname, settings, plus, minus, gff_dict, path_analysis]
        arguments.append(argument)
    
    #ribo_util.multiprocess(run_asymmetry, arguments, threads)
    print "Finished asymmetry analysis at " + str(datetime.now())

    
def run_genelist(fname, settings, plus, minus, gff, path_analysis): 
    
    '''Import settings, and assign variables'''
    
    minlength = settings['minlength']
    maxlength = settings['maxlength']
    lengthindex = range(minlength, maxlength+1)
    
    '''Load Data Files'''
    
    density_plus  = plus
    density_minus = minus 
    density_plus_rpm, density_minus_rpm = ribo_util.get_density_rpm(density_plus, density_minus, minlength, maxlength)
    
    # Parse GFF into component datasets
    
    gff_dict = gff
    alias_list       = gff_dict['Alias'] 
    strand_list      = gff_dict['Strand'] 
    start_list       = gff_dict['Start'] 
    start_codon_list = gff_dict['Start_Codon'] 
    stop_list        = gff_dict['Stop'] 
    stop_codon_list  = gff_dict['Stop_Codon'] 
    #annotation_list  = gff_dict['Annotation']
    #SD_affinity      = gff_dict['SD_affinity'] 
    
    '''Define genelist datastructure'''
    
    genelist    = {}
    rpkm_list   = []
    length_list = []
        
    '''For genes in genelist, get RPKM'''
    
    for alias, start, stop, strand in itertools.izip(alias_list, start_list,stop_list, strand_list):  
        
        # define gene start and stop 
        # RPKM calculated by excluding start peaks or stop peaks
        
        if strand == '+':
            startp = start + 30
            stopp  = stop  - 12
            genelength = stop - start
            
        elif strand == '-':
            startp = start - 30
            stopp  = stop  + 12
            genelength = start - stop
            
        # get rpkm value
        
        rpkm = ribo_util.get_RPKM(alias, startp, stopp, strand, density_plus_rpm, density_minus_rpm, minlength, maxlength)
        
        # append to list
        
        rpkm_list.append(rpkm)
        length_list.append(genelength)
        
    '''Add information to genelist dictionary and save as pickle and .csv'''
    
    genelist['Alias']       = alias_list
    #genelist['SD_Affinity'] = SD_affinity
    genelist['Start_Codon'] = start_codon_list
    genelist['Stop_Codon']  = stop_codon_list
    genelist['RPKM']        = rpkm_list
    genelist['Genelength']  = length_list
    #genelist['Annotation']  = annotation_list

    
    df = pd.DataFrame(genelist)
    df.to_csv(path_analysis + fname + '_genelist.csv')
    
    ribo_util.makePickle(genelist, path_analysis + 'genelist')

    
def genelist(inputs, paths_out, settings, gff_dict, plus_dict, minus_dict):
    files     = inputs['files']
    threads   = inputs['threads'] 
    arguments = []
    
    if not files:
        print("There are no files")
        return
    
    print "Started genelist analysis at " + str(datetime.now())

    for fname in files:
        path_analysis = paths_out['path_analysis'] + fname + '/'
        if not os.path.exists(path_analysis):
            os.makedirs(path_analysis)
        plus  = plus_dict[fname]
        minus = minus_dict[fname] 
        
        run_genelist(fname, settings, plus, minus, gff_dict, path_analysis)
        argument = [fname, settings, plus, minus, gff_dict, path_analysis]
        
        arguments.append(argument)
    
    #ribo_util.multiprocess(run_asymmetry, arguments, threads)
    print "Finished genelist analysis at " + str(datetime.now())

def run_TE(ribo_fname, seq_fname, paths_out):
    
    '''Import settings, and assign variables'''
    
    path_analysis = paths_out['path_analysis'] 
    Ribo_seq = ribo_fname
    RNA_seq  = seq_fname
    
    '''Load Data Files'''
    
    Ribo_genelist_path = path_analysis + Ribo_seq + '/' + 'genelist'
    RNA_genelist_path  = path_analysis + RNA_seq + '/' + 'genelist'
    
    Ribo_genelist = ribo_util.unPickle(Ribo_genelist_path) 
    RNA_genelist  = ribo_util.unPickle(RNA_genelist_path)
    
    if not os.path.exists(Ribo_genelist_path):
        print 'First make genelist file for ' + Ribo_seq
        return
    if not os.path.exists(Ribo_genelist_path):
        print 'First make genelist file for ' + RNA_seq
        return
    
    
    # Parse alias and RPKM info from genelists:
    
    
    Ribo_alias_list = Ribo_genelist['Alias']  
    RNA_alias_list  = RNA_genelist['Alias']
    Ribo_RPKM_list = Ribo_genelist['RPKM']  
    RNA_RPKM_list  = RNA_genelist['RPKM']  
    
    '''Define TE_list, to be appended to genelist'''
    
    TE_list = []
    
    '''Iterate through Ribo-seq alias and RPKM, 
    calculate TE for shared genes'''
    
    #first check if alias is in RNA-seq genelist:
    
    for Ribo_alias, Ribo_RPKM in itertools.izip(Ribo_alias_list, Ribo_RPKM_list):
        if Ribo_alias in RNA_alias_list:
            
            # if alias is shared, get RNA-seq RPKM value:
            
            alias_index = RNA_alias_list.index(Ribo_alias)
            RNA_RPKM    = RNA_RPKM_list[alias_index]
            
        else:
            
            # if alias not shared, set RNA_RPKM as 0 
            
            RNA_RPKM = 0
        
        # calculate TE, excluding genes with 0 RPKM for RNAseq
        
        if RNA_RPKM > 0:
            TE = Ribo_RPKM / RNA_RPKM
        else:
            TE = 0 
            
        # append TE value to list:
        
        TE_list.append(TE)
        
    '''Append TE_list to Ribo-seq genelist'''
    
    Ribo_genelist['TE'] = TE_list
            
    df = pd.DataFrame(Ribo_genelist)
    df.to_csv(path_analysis + Ribo_seq + '/' + 'genelist.csv')
    
    ribo_util.makePickle(Ribo_genelist, Ribo_genelist_path)
    
    print "Added TE value to " + Ribo_seq


    return


def run_read_composition_context(fname, inputs, settings, plus, minus, gff, path_out):
    
    minlength    = settings['minlength']
    maxlength    = settings['maxlength']
    gff_extra    = settings['gff_extra']   # extra nucleotides in gff_dict sequence (UTR sequence)
    alignment    = settings['alignment']
    windowlength = settings['windowlength']
    
    start_trim = settings['start_trim'] 
    stop_trim  = settings['stop_trim']
    
    
    density_plus  = plus
    density_minus = minus 
    gff_dict = gff
    
    alias_list  = gff_dict['Alias'] 
    strand_list = gff_dict['Strand'] 
    start_list  = gff_dict['Start'] 
    stop_list   = gff_dict['Stop']
    seq_list    = gff_dict['Sequence']
    
    positionindex = range(0, windowlength+1) 
    lengthindex   = range(minlength, maxlength+1)
    length_data   = {length : 0  for length in lengthindex}

    G_data  = {length : [0]*(windowlength+1) for length in lengthindex}
    C_data  = {length : [0]*(windowlength+1) for length in lengthindex}
    A_data  = {length : [0]*(windowlength+1) for length in lengthindex}
    U_data  = {length : [0]*(windowlength+1) for length in lengthindex}
        
    for alias, start, stop, strand, sequence in itertools.izip(alias_list, start_list,
                                                                stop_list, strand_list, seq_list):
            
        if strand == '+':
            density_dict  = density_plus
            density_start = start + start_trim
            density_stop  = stop - stop_trim          
            period = 1
            
        elif strand == '-':
            density_dict  = density_minus
            density_start = start - start_trim 
            density_stop  = stop + stop_trim
            period = -1
        
        
        #make empty density dict for the gene
        genelength    = abs(density_start - density_stop) + 1
        gene_density  = {}
        total_density = [0] * genelength
        
        #if gene is shorter than window of analysis:
        if genelength < windowlength:
            continue
        
        # gives density encompassed by seq, plus extra defined by plotlength
        for length in lengthindex:
            length_density       = density_dict[length][density_start: density_stop: period]
            gene_density[length] = length_density
            total_density        = [x + y for x, y in itertools.izip(total_density, length_density)]
        
        total_density = sum(total_density)
        
        total_read_count = 0 
                                                      
        for length in lengthindex:
            
            generated_read_list = []
            position = 0 
            
            for read_counts in gene_density[length]:
                
                
                position_3 = position + gff_extra + start_trim + 5
                position_5 = position_3 - windowlength
                newread = sequence[position_5: position_3]
                newread_list = [newread] * read_counts
                
                generated_read_list += newread_list
                
                length_data[length] += read_counts 
                total_read_count += read_counts
                
                position += 1
        
            if not generated_read_list:
                continue
                
            for generated_read in generated_read_list:
                
                read_length = len(generated_read) 
                if read_length < windowlength:
                    print 'small'
                    
                for position2 in range(1, read_length+1):
                    if generated_read[-position2] == 'G':
                        G_data[length][position2-1] += 1
                    elif generated_read[-position2] == 'C':
                        C_data[length][position2-1] += 1
                    elif generated_read[-position2] == 'A':
                        A_data[length][position2-1] += 1
                    elif generated_read[-position2] == 'T':
                        U_data[length][position2-1] += 1
                                
    print total_read_count, length_data
        
    for length in lengthindex:  
        if length_data[length] >=1:
            T = length_data[length]
        else:
            T = 1
            
        for position in range(0, windowlength+1):
                
            g = float(G_data[length][position]) / float(T) *100
            c = float(C_data[length][position]) / float(T) *100 
            a = float(A_data[length][position]) / float(T) *100 
            u = float(U_data[length][position]) / float(T) *100 
            
            G_data[length][position] = g
            C_data[length][position] = c
            A_data[length][position] = a
            U_data[length][position] = u
    
    sum_reads = sum(length_data.values())
    sum_reads = float(sum_reads)
    
    '''# calculate read length distribution: 
    
    for length in length_dist.keys():
        sum_lenght = float(length_data[length])
        read_fraction = sum_lenght / float(sum_reads) * 100
        length_dist[length] = read_fraction'''
        
    ribo_util.makePickle(G_data, path_out + 'comp_genome_G')
    ribo_util.makePickle(C_data, path_out + 'comp_genome_C')
    ribo_util.makePickle(A_data, path_out + 'comp_genome_A')
    ribo_util.makePickle(U_data, path_out + 'comp_genome_U')       
    #ribo_util.makePickle(length_dist, path_out + 'read_distribution_gene')
    
    #print read_count, pos_count, index_count
    #print excluded_genes

                  
        
def read_context(inputs, paths_out, settings, gff_dict, plus_dict, minus_dict):
    files     = inputs['files']
    threads   = inputs['threads'] 
    arguments = []
    
    if not files:
        print("There are no files")
        return
    
    print "Started  at " + str(datetime.now())

    for fname in files:
        path_analysis = paths_out['path_analysis'] + fname + '/reads_position/'
        if not os.path.exists(path_analysis):
            os.makedirs(path_analysis)
        plus  = plus_dict[fname]
        minus = minus_dict[fname] 
        
        run_read_composition_context(fname, inputs, settings, plus, minus, gff_dict, path_analysis)
        argument = [fname, inputs, plus, minus, gff_dict, path_analysis]
        arguments.append(argument)
    
    #ribo_util.multiprocess(run_read_extended, arguments, threads)
    print "Finished  at " + str(datetime.now())

    

    
def make_motif_list(motif_length):
    
    '''make list of all possible nucleotide motifs'''
        
    bases  = ['A','T','G','C']
    length = motif_length 
    
    motifs = [''.join(p) for p in itertools.product(bases, repeat = length)]
    
    return motifs

def run_motif_pausescore(fname, settings, plus, minus, gff, path_pausescore):
    
    '''define variables'''
    
    motif_length = settings['motif_length']
    minlength    = settings['minlength']
    maxlength    = settings['maxlength']
    lengthindex  = range(minlength, maxlength + 1)
    
    A_site = settings['A_site shift']
    P_site = A_site - 3
    E_site = A_site - 6
    
    frameshift      = settings['frameshift']

    plot_upstream   = settings['plot_upstream'] / 3 * 3        #change window to interval of 3
    plot_downstream = settings['plot_downstream'] / 3 * 3
    
    # define plot length
    plotlength = plot_upstream + plot_downstream + 1
    positionindex = range(0, plotlength)
    
    # load density files
    density_plus  = plus
    density_minus = minus 
    
    # load annotation
    gff_dict = gff
    
    alias_list  = gff_dict['Alias'] 
    strand_list = gff_dict['Strand'] 
    start_list  = gff_dict['Start'] 
    stop_list   = gff_dict['Stop'] 
    seq_list    = gff_dict['Sequence'] 
 
    gff_extra  = settings['gff_extra']   # extra nucleotides in gff_dict sequence (UTR sequence)
    start_trim = settings['start_trim'] / 3 * 3
    stop_trim  = settings['stop_trim'] / 3 * 3   
    
    
    minlength_1       = str(minlength)      +'_'
    maxlength_1       = str(maxlength)      +'_'
    plot_upstream_1   = str(plot_upstream)  +'_'
    plot_downstream_1 = str(plot_downstream)+'_'
    start_trim_1      = str(start_trim)     +'_'
    stop_trim_1       = str(stop_trim)      +'_'
    frameshift_1      = str(frameshift)     +'_'
    motif_length_1    = str(motif_length)   +'_'
    
    name_settings =  motif_length_1+plot_upstream_1+plot_downstream_1
    name_settings += start_trim_1+stop_trim_1+minlength_1+maxlength_1+frameshift_1
    
    # import genetic code
    aa_code, codon_code = ribo_util.get_genetic_code()
    
    '''output data structure'''
    
    # create empty datastructures to store density info:
    
    codon_motif_avgplot    = {}
    codon_motif_count      = {}
    codon_motif_score      = {}
    codon_motif_score_df   = {}
    codon_motif_list       = []
    codon_motif_count_list = []
    codon_motif_A_list     = []
    codon_motif_P_list     = []
    codon_motif_E_list     = []
    
    aa_motif_avgplot    = {}
    aa_motif_count      = {}
    aa_motif_score      = {}
    aa_motif_score_df   = {}
    aa_motif_list       = []
    aa_motif_count_list = []
    aa_motif_A_list     = []
    aa_motif_P_list     = []
    aa_motif_E_list     = []
    
    
    # generate all possible motifs of a specific length:
    motifs = make_motif_list(motif_length)
    
    for motif in motifs:
        
        # translate into an aa_motif
        aa_motif = ''
        for codon in textwrap.wrap(motif, 3):
            aa = codon_code[codon]
            aa_motif += aa
            
            '''#check for stop codons
            if aa == '_':
                stop_codon = True
            else:
                stop_codon = False
        
        # Ignore stop codon motifs
        if stop_codon == True:
            continue'''
        
        
        codon_motif_avgplot[motif] = {length : [0]*(plotlength) for length in lengthindex}
        codon_motif_count[motif]   = 0
        codon_motif_score[motif]   = [0, 0, 0] # [Asite, P site, E site] 
        
        aa_motif_avgplot[aa_motif] = {length : [0]*(plotlength) for length in lengthindex}
        aa_motif_count[aa_motif]   = 0
        aa_motif_score[aa_motif]   = [0, 0, 0] # [Asite, P site, E site] 
        
    
    '''genes in data''' 
    
    # count genes excluded from data  = [count, [names of genes]]   
    excluded_genes = {}
    excluded_genes['short']         = [0, []]
    excluded_genes['low_density']   = [0, []]
    excluded_genes['not_divisible'] = [0, []]
    # count included genes
    included_genes = [0, []]
        
        
    '''iterate through every annotated gene to get codon density info:'''
    
    for alias, start, stop, strand, sequence in itertools.izip(alias_list, start_list,stop_list, strand_list, seq_list):  
                
        ''' define start and stop positions for codons to analyze:
        
        # = codon 
        #################################### = GFF sequence  (50 extra nt)
           ##############################    = AUG to UGA
             ##########################      = density to analyze : remove start and stop peaks
                 ##################          = codons to analyze : remove plot window
        
        '''
        
        # First, define density without start and stop peaks:
        if strand == '+':
            density_dict  = density_plus
            density_start = start + start_trim + frameshift
            density_stop  = stop  - stop_trim + frameshift
            
            period = 1
            
        elif strand == '-':
            density_dict  = density_minus
            density_start = start - start_trim - frameshift
            density_stop  = stop  + stop_trim - 3 + frameshift
            
            period = -1
        
        # GFF seq has 50 extra nucleotides, so remove:
        # Also remove several codons from start and stop positions:
        
        seq_start = gff_extra + start_trim + plot_upstream + frameshift
        seq_stop  = -gff_extra - stop_trim - plot_downstream + frameshift
            
        seq       = sequence[seq_start : seq_stop]
        seqlength = len(seq)
    
        # make a list of codons in the sequence:
        motifs_in_seq = [seq[i:i+motif_length] for i in range(0, seqlength, 3)]
        
        #make empty density dict for the gene
        genelength    = abs(density_start - density_stop) + 1
        gene_density  = {}
        total_density = [0] * genelength
        
        # fill density dict with density info
        # gives density encompassed by seq, plus extra defined by plotlength
        for length in lengthindex:
            length_density       = density_dict[length][density_start: density_stop: period]
            length_density_float = [float(i) for i in length_density]
            gene_density[length] = length_density
            total_density        = [x + y for x, y in itertools.izip(total_density, length_density)]
                        
        gene_avgreads = float(sum(total_density)) / float(genelength)
        # exclude genes that are not divisable by 3
        if seqlength % 3 != 0:
            excluded_genes['not_divisible'][0] += 1
            excluded_genes['not_divisible'][1].append(alias)
            continue
        
        # exclude genes shorter than plot
        if seqlength < plotlength + 1:
            excluded_genes['short'][0] += 1
            excluded_genes['short'][1].append(alias)
            continue
            
        elif seq_start - seq_stop > len(sequence):
            excluded_genes['short'][0] += 1
            excluded_genes['short'][1].append(alias)
            continue
            
        # normalize density by total gene density
        if gene_avgreads < 0.1:
            excluded_genes['low_density'][0] += 1
            excluded_genes['low_density'][1].append(alias)
            continue
            
        else: 
            relative_density = {}
            
            for length in lengthindex:
                relative_density[length] = [reads / gene_avgreads for reads in gene_density[length]]
        
        # add data to dataframe
        motif_index = 0
        for motif in motifs_in_seq:
            if len(motif) < motif_length:
                continue
            motif_position = (motif_index * 3)
            
            aa_motif = ''
            for codon in textwrap.wrap(motif, 3):
                aa = codon_code[codon]
                aa_motif += aa
            
            motif_index        += 1
            codon_motif_count[motif] += 1
            aa_motif_count[aa_motif] += 1
            
            for length in lengthindex:
                for position in positionindex:
                    
                    density_position = motif_position + position
                    density          = relative_density[length][density_position]
                    
                    codon_motif_avgplot[motif][length][position] += density
                    aa_motif_avgplot[aa_motif][length][position] += density
                    
        included_genes[0] += 1
        included_genes[1].append(alias)
    
    #divide data by total instances of the codon_motif or aa_motif
    for codon_motif in codon_motif_avgplot.keys():
        for length in lengthindex:
            for position in positionindex:
                if codon_motif_count[codon_motif] == 0:
                    continue 
                else: 
                    codon_motif_avgplot[codon_motif][length][position] = codon_motif_avgplot[codon_motif][length][position] / codon_motif_count[codon_motif]
    
    for aa_motif in aa_motif_avgplot.keys():
        for length in lengthindex:
            for position in positionindex:
                if aa_motif_count[aa_motif] == 0:
                    continue
                else: 
                    aa_motif_avgplot[aa_motif][length][position] = aa_motif_avgplot[aa_motif][length][position] / aa_motif_count[aa_motif]
    
    '''Convert data for plotting and csv'''
    
    codon_motif_data     = {}
    codon_motif_data_sum = {}
    aa_motif_data     = {}
    aa_motif_data_sum = {}
    
    A_site_shift = plot_upstream - A_site 
    P_site_shift = plot_upstream - P_site 
    E_site_shift = plot_upstream - E_site 
    
    for codon_motif in codon_motif_avgplot.keys():
        df = pd.DataFrame(codon_motif_avgplot[codon_motif])
        df = df.T
        df = df.reindex(index=df.index[::-1])

        plot_range  = len(df.columns)
        shift_start = - plot_upstream
        shift_stop  = plot_range - plot_upstream

        df.columns = pd.RangeIndex(start = shift_start, stop = shift_stop, step = 1)

        codon_motif_data[codon_motif]     = df
        codon_motif_data_sum[codon_motif] = df.sum(0)
        
        codon_motif_score[codon_motif][0] = sum(codon_motif_data_sum[codon_motif][A_site_shift: A_site_shift + 3:1])/3
        codon_motif_score[codon_motif][1] = sum(codon_motif_data_sum[codon_motif][P_site_shift: P_site_shift + 3:1])/3 
        codon_motif_score[codon_motif][2] = sum(codon_motif_data_sum[codon_motif][E_site_shift: E_site_shift + 3:1])/3
        
        count = codon_motif_count[codon_motif]
        
        codon_motif_list.append(codon_motif)
        codon_motif_count_list.append(count)
        codon_motif_A_list.append(codon_motif_score[codon_motif][0])
        codon_motif_P_list.append(codon_motif_score[codon_motif][1])
        codon_motif_E_list.append(codon_motif_score[codon_motif][2])
        
    for aa_motif in aa_motif_avgplot.keys():
        df = pd.DataFrame(aa_motif_avgplot[aa_motif])
        df = df.T
        df = df.reindex(index=df.index[::-1])

        plot_range  = len(df.columns)
        shift_start = - plot_upstream
        shift_stop  = plot_range - plot_upstream

        df.columns = pd.RangeIndex(start = shift_start, stop = shift_stop, step = 1)

        aa_motif_data[aa_motif]     = df
        aa_motif_data_sum[aa_motif] = df.sum(0)
        
        aa_motif_score[aa_motif][0] = sum(aa_motif_data_sum[aa_motif][A_site_shift: A_site_shift + 3:1])/3
        aa_motif_score[aa_motif][1] = sum(aa_motif_data_sum[aa_motif][P_site_shift: P_site_shift + 3:1])/3                
        aa_motif_score[aa_motif][2] = sum(aa_motif_data_sum[aa_motif][E_site_shift: E_site_shift + 3:1])/3
        
        count = aa_motif_count[aa_motif]
        
        aa_motif_list.append(aa_motif)
        aa_motif_count_list.append(count)
        aa_motif_A_list.append(aa_motif_score[aa_motif][0])
        aa_motif_P_list.append(aa_motif_score[aa_motif][1])
        aa_motif_E_list.append(aa_motif_score[aa_motif][2])
    
    codon_motif_score_df['Motif']  = codon_motif_list
    codon_motif_score_df['A_site'] = codon_motif_A_list
    codon_motif_score_df['P_site'] = codon_motif_P_list
    codon_motif_score_df['E_site'] = codon_motif_E_list
    codon_motif_score_df['Count']  = codon_motif_count_list
    
    aa_motif_score_df['Motif']  = aa_motif_list
    aa_motif_score_df['A_site'] = aa_motif_A_list
    aa_motif_score_df['P_site'] = aa_motif_P_list
    aa_motif_score_df['E_site'] = aa_motif_E_list
    aa_motif_score_df['Count']  = aa_motif_count_list
    
    codon_motif_df = pd.DataFrame(codon_motif_score_df)
    codon_motif_df.to_csv(path_pausescore + 'codon_motif_scores_' + name_settings + '.csv')
    
    aa_motif_df = pd.DataFrame(aa_motif_score_df)
    aa_motif_df.to_csv(path_pausescore + 'aa_motif_scores_' + name_settings + '.csv')

    
    '''ribo_util.makePickle(motif_data,  path_pausescore + name_settings +'motif_HM_data' protocol=pickle.HIGHEST_PROTOCOL) 
    ribo_util.makePickle(motif_data_sum, path_pausescore + name_settings +'motif_plot_data', protocol=pickle.HIGHEST_PROTOCOL)  
    ribo_util.makePickle(motif_score_df, path_pausescore + name_settings +'motif_scores'   , protocol=pickle.HIGHEST_PROTOCOL)
    '''
    
    print excluded_genes['short'][0]    
    print excluded_genes['low_density'][0]
    print excluded_genes['not_divisible'][0]
    print included_genes[0]
    return 


def motif_pausescore(inputs, paths_out, settings, gff_dict, plus_dict, minus_dict):
    
    files     = inputs['files']
    threads   = inputs['threads'] 
    multi     = 'no'
    arguments = []
    
    if not files:
        print("There are no files")
        return
    
    print "Started pause score analysis at " + str(datetime.now())

    for fname in files:
        path_pausescore = paths_out['path_analysis'] + fname + '/pause_score/'
        if not os.path.exists(path_pausescore):
            os.makedirs(path_pausescore)
        plus  = plus_dict[fname]
        minus = minus_dict[fname] 
       
        if not multi == 'yes':
            run_motif_pausescore(fname, settings, plus, minus, gff_dict, path_pausescore)
        else:     
            argument = [fname, settings, plus, minus, gff_dict, path_pausescore]
            arguments.append(argument)
    
    if multi == 'yes':
        ribo_util.multiprocess(run_motif_pausescore, arguments, threads)
    
    print "Finished pause score analysis at " + str(datetime.now())
    
    return

def run_pausescore_waves(fname, settings, plus, minus, gff, path_pausescore):
    
    '''define variables'''
    
    minlength   = settings['minlength']
    maxlength   = settings['maxlength']
    lengthindex = range(minlength, maxlength + 1)
    
    A_site = settings['A_site shift']
    P_site = A_site - 3
    E_site = A_site - 6
    
    frameshift      = settings['frameshift']

    plot_upstream   = settings['plot_upstream_wave'] / 3 * 3        #change window to interval of 3
    plot_downstream = settings['plot_downstream_wave'] / 3 * 3
    
    next_codon = settings['next_codon']
    
    # define plot length
    plotlength = plot_upstream + plot_downstream + 1
    positionindex = range(0, plotlength)
    
    # load density files
    density_plus  = plus
    density_minus = minus 
    
    # load annotation
    gff_dict = gff
    
    alias_list  = gff_dict['Alias'] 
    strand_list = gff_dict['Strand'] 
    start_list  = gff_dict['Start'] 
    stop_list   = gff_dict['Stop'] 
    seq_list    = gff_dict['Sequence'] 
 
    gff_extra  = settings['gff_extra']   # extra nucleotides in gff_dict sequence (UTR sequence)
    start_trim = settings['start_trim'] / 3 * 3
    stop_trim  = settings['stop_trim'] / 3 * 3   
    
    
    minlength_1       = str(minlength)      +'_'
    maxlength_1       = str(maxlength)      +'_'
    plot_upstream_1   = str(plot_upstream)  +'_'
    plot_downstream_1 = str(plot_downstream)+'_'
    start_trim_1      = str(start_trim)     +'_'
    stop_trim_1       = str(stop_trim)      +'_'
    frameshift_1      = str(frameshift)     +'_'
    next_codon_1      =  str(next_codon)    +'_'
    
    name_settings = plot_upstream_1+plot_downstream_1+start_trim_1+stop_trim_1
    name_settings += minlength_1+maxlength_1+frameshift_1+next_codon_1
    
    
    
    # import genetic code
    aa_code, codon_code = ribo_util.get_genetic_code()
    
    
    '''output data structure'''
    
    # aa/codon avgplot      = { aa: {length: [values]}} 
    # aa/codon count        = { aa: N}
    # aa/codon score        = { aa: [aa_list], A site score: [A_list]...
                                
    aa_avgplot  = {}
    aa_count    = {}
    aa_score    = {}
    aa_score_df = {}
    aa_list     = []
    aa_A_list   = []
    aa_P_list   = []
    aa_E_list   = []
    
    codon_avgplot  = {}
    codon_count    = {}
    codon_score    = {}
    codon_score_df = {}
    codon_list     = []
    codon_A_list   = []
    codon_P_list   = []
    codon_E_list   = []
    
    # create empty datastructures to store density info:
    for aa in aa_code.keys():
        aa_avgplot[aa] = {length : [0]*(plotlength) for length in lengthindex}
        aa_count[aa]   = 0
        aa_score[aa]   = [0, 0, 0] # [Asite, P site, E site] 
        
    for codon in codon_code.keys():
        codon_avgplot[codon] = {length : [0]*(plotlength) for length in lengthindex}
        codon_count[codon]   = 0
        codon_score[codon]   = [0, 0, 0] # [Asite, P site, E site] 
    
    '''genes in data''' 
    
    # count genes excluded from data  = [count, [names of genes]]   
    excluded_genes = {}
    excluded_genes['short']         = [0, []]
    excluded_genes['low_density']   = [0, []]
    excluded_genes['not_divisible'] = [0, []]
    # count included genes
    included_genes = [0, []]
        
        
    '''iterate through every annotated gene to get codon density info:'''
    
    for alias, start, stop, strand, sequence in itertools.izip(alias_list, start_list,stop_list, strand_list, seq_list):  
                
        ''' define start and stop positions for codons to analyze:
        
        # = codon 
        #################################### = GFF sequence  (50 extra nt)
           ##############################    = AUG to UGA
             ##########################      = density to analyze : remove start and stop peaks
                 ##################          = codons to analyze : remove plot window
        
        '''
        
        # First, define density without start and stop peaks:
        if strand == '+':
            density_dict  = density_plus
            density_start = start + start_trim + frameshift
            density_stop  = stop  - stop_trim + frameshift
            
            period = 1
            
        elif strand == '-':
            density_dict  = density_minus
            density_start = start - start_trim - frameshift
            density_stop  = stop  + stop_trim - 3 + frameshift
            
            period = -1
        
        # GFF seq has 50 extra nucleotides, so remove:
        # Also remove several codons from start and stop positions:
        
        codon_seq_start = gff_extra + start_trim + plot_upstream + frameshift
        codon_seq_stop  = -gff_extra - stop_trim - plot_downstream + frameshift
            
        seq       = sequence[codon_seq_start : codon_seq_stop]
        seqlength = len(seq)
    
        #make empty density dict for the gene
        genelength    = abs(density_start - density_stop) + 1
        gene_density  = {}
        total_density = [0] * genelength
        
        # fill density dict with density info
        # gives density encompassed by seq, plus extra defined by plotlength
        for length in lengthindex:
            length_density       = density_dict[length][density_start: density_stop: period]
            length_density_float = [float(i) for i in length_density]
            gene_density[length] = length_density
            total_density        = [x + y for x, y in itertools.izip(total_density, length_density)]
            
        gene_avgreads = float(sum(total_density)) / float(genelength)
        # exclude genes that are not divisable by 3
        if seqlength % 3 != 0:
            excluded_genes['not_divisible'][0] += 1
            excluded_genes['not_divisible'][1].append(alias)
            continue
        
        # exclude genes shorter than plot
        if seqlength < plotlength + 1:
            excluded_genes['short'][0] += 1
            excluded_genes['short'][1].append(alias)
            continue
            
        elif codon_seq_start - codon_seq_stop > len(sequence):
            excluded_genes['short'][0] += 1
            excluded_genes['short'][1].append(alias)
            continue
            
        # normalize density by total gene density
        if gene_avgreads < 0.1:
            excluded_genes['low_density'][0] += 1
            excluded_genes['low_density'][1].append(alias)
            continue
            
        else: 
            relative_density = {}
            
            for length in lengthindex:
                relative_density[length] = [reads / gene_avgreads for reads in gene_density[length]]
        
        # add data to dataframe
        codon_index = 0
        
        # make a list of codons in the sequence:
        codons_seq = [seq[i:i+3] for i in range(0, seqlength, 3)]
        aa_seq     = [codon_code[codon] for codon in codons_seq]
        
        for codon in codons_seq:
            codon_position = (codon_index * 3)
            aa = aa_seq[codon_index]
            
            if next_codon == 'yes':
                if aa in aa_seq[codon_index+1:plot_downstream + 1]:
                    codon_index += 1
                    continue
            
            codon_index        += 1
            codon_count[codon] += 1
            aa_count[aa]       += 1
            
            for length in lengthindex:
                for position in positionindex:
                    
                    density_position = codon_position + position
                    density          = relative_density[length][density_position]
                    
                    codon_avgplot[codon][length][position] += density
                    aa_avgplot[aa][length][position]       += density
                    
        included_genes[0] += 1
        included_genes[1].append(alias)
    
    #divide data by total instances of the codon or aa
    for codon in codon_avgplot.keys():
        for length in lengthindex:
            for position in positionindex:
                if codon_count[codon] == 0:
                    continue 
                else: 
                    codon_avgplot[codon][length][position] = codon_avgplot[codon][length][position] / codon_count[codon]
    
    for aa in aa_avgplot.keys():
        for length in lengthindex:
            for position in positionindex:
                if aa_count[aa] == 0:
                    continue
                else: 
                    aa_avgplot[aa][length][position] = aa_avgplot[aa][length][position] / aa_count[aa]
                
    '''Convert data for plotting and csv'''
    
    codon_data     = {}
    codon_data_sum = {}
    aa_data     = {}
    aa_data_sum = {}
    
    A_site_shift = plot_upstream - A_site 
    P_site_shift = plot_upstream - P_site 
    E_site_shift = plot_upstream - E_site 
    
    for codon in codon_avgplot.keys():
        df = pd.DataFrame(codon_avgplot[codon])
        df = df.T
        df = df.reindex(index=df.index[::-1])

        plot_range  = len(df.columns)
        shift_start = - plot_upstream
        shift_stop  = plot_range - plot_upstream

        df.columns = pd.RangeIndex(start = shift_start, stop = shift_stop, step = 1)

        codon_data[codon]     = df
        codon_data_sum[codon] = df.sum(0)
        
        codon_score[codon][0] = sum(codon_data_sum[codon][A_site_shift: A_site_shift + 3:1])/3
        codon_score[codon][1] = sum(codon_data_sum[codon][P_site_shift: P_site_shift + 3:1])/3                            
        codon_score[codon][2] = sum(codon_data_sum[codon][E_site_shift: E_site_shift + 3:1])/3
        
        codon_list.append(codon)
        codon_A_list.append(codon_score[codon][0])
        codon_P_list.append(codon_score[codon][1])
        codon_E_list.append(codon_score[codon][2])
        
    for aa in aa_avgplot.keys():
        df = pd.DataFrame(aa_avgplot[aa])
        df = df.T
        df = df.reindex(index=df.index[::-1])

        plot_range  = len(df.columns)
        shift_start = - plot_upstream
        shift_stop  = plot_range - plot_upstream

        df.columns = pd.RangeIndex(start = shift_start, stop = shift_stop, step = 1)

        aa_data[aa]     = df
        aa_data_sum[aa] = df.sum(0) 
                                    
        aa_score[aa][0] = sum(aa_data_sum[aa][A_site_shift: A_site_shift + 3:1])/3    
        aa_score[aa][1] = sum(aa_data_sum[aa][P_site_shift: P_site_shift + 3:1])/3
        aa_score[aa][2] = sum(aa_data_sum[aa][E_site_shift: E_site_shift + 3:1])/3  
        
        aa_list.append(aa)
        aa_A_list.append(aa_score[aa][0])
        aa_P_list.append(aa_score[aa][1])
        aa_E_list.append(aa_score[aa][2])
    
    codon_score_df['Codon']  = codon_list
    codon_score_df['A_site'] = codon_A_list
    codon_score_df['P_site'] = codon_P_list
    codon_score_df['E_site'] = codon_E_list
    
    aa_score_df['Amino Acid'] = aa_list
    aa_score_df['A_site']     = aa_A_list
    aa_score_df['P_site']     = aa_P_list
    aa_score_df['E_site']     = aa_E_list
    
    
    codon_df = pd.DataFrame(codon_score_df)
    codon_df.to_csv(path_pausescore + 'codon_scores.csv')
    
    aa_df = pd.DataFrame(aa_score_df)
    aa_df.to_csv(path_pausescore + 'aa_scores.csv')
    
    ribo_util.makePickle(codon_data,     path_pausescore + name_settings +'codon_HM_data'  , protocol=pickle.HIGHEST_PROTOCOL) 
    ribo_util.makePickle(codon_data_sum, path_pausescore + name_settings +'codon_plot_data', protocol=pickle.HIGHEST_PROTOCOL)  
    ribo_util.makePickle(codon_score_df, path_pausescore + name_settings +'codon_scores'   , protocol=pickle.HIGHEST_PROTOCOL)
    ribo_util.makePickle(aa_data,        path_pausescore + name_settings +'aa_HM_data'     , protocol=pickle.HIGHEST_PROTOCOL) 
    ribo_util.makePickle(aa_data_sum,    path_pausescore + name_settings +'aa_plot_data'   , protocol=pickle.HIGHEST_PROTOCOL)
    ribo_util.makePickle(aa_score_df,    path_pausescore + name_settings +'aa_scores'      , protocol=pickle.HIGHEST_PROTOCOL)
    
    print excluded_genes['short'][0]    
    print excluded_genes['low_density'][0]
    print excluded_genes['not_divisible'][0]
    print included_genes[0]
    return 


def pausescore_waves(inputs, paths_out, settings, gff_dict, plus_dict, minus_dict):
    
    files     = inputs['files']
    threads   = inputs['threads'] 
    multi     = inputs['multiprocess']
    arguments = []
    
    if not files:
        print("There are no files")
        return
    
    print "Started pause score analysis at " + str(datetime.now())

    for fname in files:
        path_pausescore = paths_out['path_analysis'] + fname + '/pause_score/waves/'
        if not os.path.exists(path_pausescore):
            os.makedirs(path_pausescore)
        plus  = plus_dict[fname]
        minus = minus_dict[fname] 
       
        if not multi == 'yes':
            run_pausescore_waves(fname, settings, plus, minus, gff_dict, path_pausescore)
        else:     
            argument = [fname, settings, plus, minus, gff_dict, path_pausescore]
            arguments.append(argument)
    
    if multi == 'yes':
        ribo_util.multiprocess(run_pausescore_waves, arguments, threads)
    
    print "Finished pause score analysis at " + str(datetime.now())
    
    return

def run_read_extended(fname, inputs, settings, plus, minus, gff, path_out):
    
    minlength = settings['minlength']
    maxlength = settings['maxlength']
    
    density_plus  = plus
    density_minus = minus 
    gff_dict = gff
    
    alias_list  = gff_dict['Alias'] 
    strand_list = gff_dict['Strand'] 
    start_list  = gff_dict['Start'] 
    stop_list   = gff_dict['Stop']
    seq_list    = gff_dict['Sequence']
    
    positionindex = range(0, maxlength+1) 
    lengthindex   = range(minlength, maxlength+1)
        
    length_data = {length : 0  for length in lengthindex}

    G_data  = {length : [0]*(maxlength+1) for length in lengthindex}
    C_data  = {length : [0]*(maxlength+1) for length in lengthindex}
    A_data  = {length : [0]*(maxlength+1) for length in lengthindex}
    U_data  = {length : [0]*(maxlength+1) for length in lengthindex}
       
    for alias, start, stop, strand, sequence in itertools.izip(
        alias_list, start_list,stop_list, strand_list, seq_list):
        
        if strand == '+':
            
            genelength = len(range(start+30, stop + 1)) 
            
            for length in lengthindex:
                for position in range(0, genelength):
                    
                    density = density_plus[length][start+30 + position] 
                    if density == 0:
                        continue
                        
                    position_end = position + 40
                    position_beg = position_end - maxlength
                    newread = sequence[position_beg: position_end]
                    newread_list = [newread] * density
                    length_data[length] += density       
                    
                    for newread in newread_list:
                        for position2 in range(1, maxlength+1):
                            if newread[-position2] == 'G':
                                G_data[length][position2-1] += 1
                            elif newread[-position2] == 'C':
                                C_data[length][position2-1] += 1
                            elif newread[-position2] == 'A':
                                A_data[length][position2-1] += 1
                            elif newread[-position2] == 'T':
                                U_data[length][position2-1] += 1
                                
    for length in lengthindex:  
        if length_data[length] >=1:
            T = length_data[length]
        else:
            T = 1
            
        for position in range(0, maxlength):
                
            g = float(G_data[length][position]) / float(T) *100
            c = float(C_data[length][position]) / float(T) *100 
            a = float(A_data[length][position]) / float(T) *100 
            u = float(U_data[length][position]) / float(T) *100 
            
            G_data[length][position] = g
            C_data[length][position] = c
            A_data[length][position] = a
            U_data[length][position] = u
    
    ribo_util.makePickle(G_data, path_out + 'comp_genome_G')
    ribo_util.makePickle(C_data, path_out + 'comp_genome_C')
    ribo_util.makePickle(A_data, path_out + 'comp_genome_A')
    ribo_util.makePickle(U_data, path_out + 'comp_genome_U')           
                  
                    
def read_extended(inputs, settings, paths_out, gff_dict, plus_dict, minus_dict):
    files     = inputs['files']
    threads   = inputs['threads'] 
    arguments = []
    
    if not files:
        print("There are no files")
        return
    
    print "Started  at " + str(datetime.now())

    for fname in files:
        path_analysis = paths_out['path_analysis'] + fname + '/'
        if not os.path.exists(path_analysis):
            os.makedirs(path_analysis)
        plus  = plus_dict[fname]
        minus = minus_dict[fname] 
        
        run_read_extended(fname, inputs, settings, plus, minus, gff_dict, path_analysis)
        argument = [fname, inputs, settings, plus, minus, gff_dict, path_analysis]
        arguments.append(argument)
    
    #ribo_util.multiprocess(run_read_extended, arguments, threads)
    print "Finished  at " + str(datetime.now())

    

    
    
    
    