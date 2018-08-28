from datetime import datetime
from multiprocessing import Process, Pool
import os, time
import subprocess
import struct
import cPickle as pickle
import csv
from BCBio import GFF
from Bio import Seq
import itertools

import pandas as pd
import numpy as np
from IPython.display import display

import ribo_util
import ribo_main
import ribo_analysis
import ribo_plot


'''These functions use known transcriptional ends from B subtilis to align desnity to those ends'''

def get_transcript_end(file_in): 
    '''function uses data from :Prediction of Transcriptional Terminators in Bacillus subtilis and Related Species 
    Michiel J. L. de Hoon, Yuko Makita, Kenta Nakai, and Satoru Miyan:
    Extracts last gene in operon, and position of transcription termination hairpin relative to stop codon
    as two lists. '''
    
    localization = pd.read_csv(file_in)
    local = localization.to_dict(orient='list')
    
    aliaslist          = local['Last gene']
    dist_from_lastgene = local['Position wrt stop codon']
    
    return aliaslist, dist_from_lastgene
    
def run_avggene_transcriptend(fname, settings, plus, minus, gff, path_start, path_stop):
    next_gene      = settings['next_gene']    
    equal_weight   = settings['equal_weight']
    length_in_ORF  = settings['length_in_ORF']
    length_out_ORF = settings['length_out_ORF']
    minlength      = settings['minlength']
    maxlength      = settings['maxlength']
    
    density_plus  = plus
    density_minus = minus 
    gff_dict = gff
    
    window_length = length_in_ORF + length_out_ORF
    positionindex = range(0, window_length) 
    lengthindex   = range(minlength, maxlength+1)
    
    alias_list  = gff_dict['Alias'] 
    strand_list = gff_dict['Strand'] 
    start_list  = gff_dict['Start'] 
    stop_list   = gff_dict['Stop'] 
    
    lastgene, dist_from_lastgene = get_transcript_end('/Volumes/HDD/data/gff/BS/rho-independent_ends.csv')

    # datastructure:
    # for heatmap - each dict has read counts at each position separated by length keys
    # for avggene summary - dictionary with position as keys 
    
    averagegene_start = {length : [0]*(window_length) for length in lengthindex}
    averagegene_stop  = {length : [0]*(window_length) for length in lengthindex}
    start_all = {position : 0 for position in positionindex}
    stop_all  = {position : 0 for position in positionindex}
    
    genes_too_close    = []  
    genes_below_thresh = []

    for alias, start, stop, strand in itertools.izip(alias_list, start_list,stop_list, strand_list):  
        
        if alias not in lastgene:
            continue
        else:
            getindex = lastgene.index(alias)
            endpos   = dist_from_lastgene[getindex]
            
        
        genelength = abs(start - stop) + 1  
        
        nextgene = ribo_util.nextgene(alias, gff_dict)
        
        # define plot start and stop window for genes in + or - strand:
        
        if strand == '+':
            stop = stop + endpos
            start_window = start - length_out_ORF 
            stop_window  = stop - length_in_ORF 
            stop_window_max = stop + length_out_ORF 
            
        if strand == '-':
            stop = stop - endpos
            start_window = start + length_out_ORF 
            stop_window  = stop + length_in_ORF 
            stop_window_max = stop - length_out_ORF 

        if next_gene > 0:
            if nextgene['distance'] < next_gene:           # exclude genes that are too close 
                genes_too_close.append(alias)
                genes_too_close.append(nextgene['alias'])
                continue
        
        if alias in genes_too_close:
            continue 
            
        if genelength < length_in_ORF + 20:            #exclude genes that are too small
            continue 
        
        if equal_weight == 'y':
            gene_reads, length_reads = ribo_util.get_genecounts(start_window, stop_window, strand, density_plus, 
                                             density_minus, minlength, maxlength)
            gene_reads = sum(gene_reads)
            if gene_reads < 10:
                genes_below_thresh.append(alias)
                continue
                
        else:
            gene_reads = 1
        
        gene_reads = float(gene_reads)        
        
        genomelength = len(density_plus[density_plus.keys()[0]])
        
        if strand == '+':
            
            if not 0 <= start_window < genomelength:
                continue
            if not 0 <= stop_window_max < genomelength:
                continue
            
            density_dict = density_plus
            for length in lengthindex:
                for position in positionindex:
                    start_density = density_dict[length][start_window + position] 
                    start_density = float(start_density) / gene_reads
                    averagegene_start[length][position] += start_density
                    start_all[position] += start_density
                    
                    stop_density = density_dict[length][stop_window + position] 
                    stop_density = float(stop_density) / gene_reads
                    averagegene_stop[length][position] += stop_density
                    stop_all[position]  += stop_density

        if strand == '-':
            
            if not 0 <= start_window < genomelength:
                continue
            if not 0 <= stop_window_max < genomelength:
                continue
                
            density_dict = density_minus
            for length in lengthindex:
                for position in positionindex:
                    start_density = density_dict[length][start_window - position]
                    start_density = float(start_density) / gene_reads
                    averagegene_start[length][position] += start_density
                    start_all[position] += start_density
                    
                    stop_density = density_dict[length][stop_window - position] 
                    stop_density = float(stop_density) / gene_reads
                    averagegene_stop[length][position] += stop_density 
                    stop_all[position]  += stop_density 
    
    genes_excluded = list(set(genes_too_close))
    print '\t' + str(len(genes_excluded)) + ' genes excluded from ' + fname + str(len(genes_below_thresh))
    
    ribo_util.makePickle(start_all, path_start + '_all_end', protocol=pickle.HIGHEST_PROTOCOL)                   
    ribo_util.makePickle(averagegene_start, path_start + '_HM_end', protocol=pickle.HIGHEST_PROTOCOL)
    ribo_util.makePickle(stop_all, path_stop + '_all_end', protocol=pickle.HIGHEST_PROTOCOL)
    ribo_util.makePickle(averagegene_stop, path_stop + '_HM_end', protocol=pickle.HIGHEST_PROTOCOL)
    
def avggenes_transcriptend(inputs, paths_out, settings, gff_dict, plus_dict, minus_dict): 
    files     = inputs['files']
    threads   = inputs['threads'] 
    arguments = []
    
    if not files:
        print("There are no files")
        return
    
    print "Started avggenes at " + str(datetime.now())

    for fname in files:
        path_analysis = paths_out['path_analysis'] + fname + '/'
        if not os.path.exists(path_analysis):
            os.makedirs(path_analysis)
        path_start = path_analysis + 'avg_start'
        path_stop  = path_analysis + 'avg_stop'
        plus  = plus_dict[fname]
        minus = minus_dict[fname] 
        
        run_avggene_transcriptend(fname, settings, plus, minus, gff_dict, path_start, path_stop)
        argument = [fname, settings, plus, minus, gff_dict, path_start, path_stop]
        arguments.append(argument)
    
    #ribo_util.multiprocess(run_avggenes, arguments, threads)
    print "Finished avggenes at " + str(datetime.now())
    
    
    '''These functions search density files to call orfs. Uses start codon peak characteristics 
        and ribosome density to identify translating regions. '''
    
def find_ORF(fname, inputs, settings, plus, minus, gff, sequence, path_start, path_stop):
    
    density_plus, density_minus  =  ribo_util.merge_density_lenghts(density_plus, density_minus, minlength, maxlength)   
    
    density_plus  = plus
    density_minus = reversed(minus) 
    gff_dict = gff
    
    lenght_genome = len(density_plus)
    
    de_novo_dict = {}
    alias = []
    startpos = []
    stoppos = [] 
    known = []
    
    for position in range(0, length_genome):
        
        if position + 21 > length_genome:         # prevents end issues 
            continue
            
        windowlenght = 21
        
        high_rpclist = []
        
        '''analyze 21 nt window of density file to find regions of ribosome occupancy'''
        for window_pos in range(position, position + windowlenght):
            
            density_window = density_plus[window_pos: window_pos + 21]     # look at windows of density
            
            sum_density = sum(density_window)
            rpc         = float(sum_density) / 7     # calculate rpc
            
            if not rpc > 0.01:
                high_rpclist.append(False)    
            else:
                high_rpclist.append(True)
                
        if rpc in high_rpclist == False:     # if region has windowed rpc less than limit, continue
            continue
            
            
        

        
    
    
    
    alias_list  = gff_dict['Alias'] 
    strand_list = gff_dict['Strand'] 
    start_list  = gff_dict['Start'] 
    stop_list   = gff_dict['Stop'] 
    
    lastgene, dist_from_lastgene = get_transcript_end('/Volumes/HDD/data/gff/BS/rho-independent_ends.csv')

    # datastructure:
    # for heatmap - each dict has read counts at each position separated by length keys
    # for avggene summary - dictionary with position as keys 
    
    averagegene_start = {length : [0]*(window_length) for length in lengthindex}
    averagegene_stop  = {length : [0]*(window_length) for length in lengthindex}
    start_all = {position : 0 for position in positionindex}
    stop_all  = {position : 0 for position in positionindex}
    
    genes_too_close    = []  
    genes_below_thresh = []

    for alias, start, stop, strand in itertools.izip(alias_list, start_list,stop_list, strand_list):  
    
    
    
    
    
    
    
