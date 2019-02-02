'''
    Script to filter, align and 3' map ribosome footprints
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


from datetime import datetime
from multiprocessing import Process
import os, time
import subprocess
import struct
import csv
from BCBio import GFF
from Bio import Seq
import cPickle as pickle


import pandas as pd
import numpy as np
from IPython.display import display

import ribo_util


'''
Table of Contents:

    -- Filtering:
        * skewer v0.2.2 
        
    -- Aligning:
        * bowtie v0.12.7
        
    -- Density:
        * 3' alignment
    

        
'''

############################
#####    Filtering     #####
############################



def run_filter(inputs, paths_in, paths_out): # all arguments = dict
    '''
    Filter reads using skewer
    '''
    
    files        = inputs['files']
    run          = inputs['run_filtering']
    minlength    = inputs['minlength']
    maxlength    = inputs['maxlength']
    phred_cutoff = inputs['phred_cutoff']
    linker       = inputs['linker']
    threads      = inputs['threads']    # filterreads has its own threading,
    filtering    = []
    log_data     = {}
    if not files:
        print("There are no files")
        return
    
    for fname in files: 
                       
        file_in  = paths_in['path_fastq'] + fname
        file_out = paths_out['path_filter'] + fname 
        file_log = paths_out['path_log'] + fname + '_filter'
                    
        if not run == 'yes':
            if not os.path.exists(file_out+'-trimmed.fastq'):
                print "ERROR: " + fname + " has not been filtered, change run setting"
                continue
            else: 
                print fname + " has been filtered"
                continue
                
        if not os.path.exists(file_in):
            print "ERROR: " + fname + " has no FASTQ file, has been removed from analysis"
            inputs['files'].remove(fname)
            continue
        
        command_to_run = 'skewer -x %s -Q %d  -l %d -L %d -o %s --quiet -t %d %s 1>>%s 2>%s' % (
            linker, 
            phred_cutoff, 
            minlength,
            maxlength, 
            file_out, 
            threads, 
            file_in, 
            file_log,
            file_log
            )
        
        #Add filter parameters to log:
        
        log_data['settings'] = {'linker': linker, 'phred_cutoff': phred_cutoff,
                                'minlength': minlength, 'maxlength': maxlength}
        
        log_function = 'ribo_density'
        ribo_util.analysis_log(fname, log_function, log_data, paths_in, paths_out)
        
        filtering.append(command_to_run)
        
    print "-----FILTER-----"
    print '\nFiles to filter: ' + ', '.join(files)
    print "Filter parameters are: \nmin length = %s \nmax length = %s \nphred cutoff = %s " % (
                    minlength, maxlength, phred_cutoff)
    print "\n\tStarted filtering at " + str(datetime.now())
    
    ribo_util.subprocess_wf(filtering, 1)
    
    print "\tFinished filtering at " + str(datetime.now())
    print "\tCOMPLETED FILTERING"
    
    return inputs
        
        
        

############################
#####     Aligning     #####   
############################

        
def run_align(inputs, paths_in, paths_out): # all arguments = dict
    '''Bowtie align'''
    
    run      = inputs['run_bowtie']
    files    = inputs['files']
    threads  = inputs['cores']     # bowtie uses 1 core per instance
    
    
    if not files:
        print("There are no files")
        return
    
    ladder     = []
    tRNA       = []
    rRNA       = []
    chromosome = []
    
    for fname in files:
        if not run == 'yes':
            if not os.path.exists(paths_out['path_chr'] + fname + '_match.SAM'):
                print "ERROR: " + fname + " has not been aligned, change run settings"
                continue
            else: 
                print fname + " has been aligned"
                continue
                
        if not os.path.exists(paths_out['path_filter'] + fname + '-trimmed.fastq'):
            print "ERROR: " + fname + " has no filtered file, has been removed from analysis"
            inputs['files'].remove(fname)
            continue
        
        file_log = paths_out['path_log'] + fname + '_bowtie'
        
        # bowtie_1 will rewrite log
        bowtie_1 =  '%s -v 2 -y -m 1 -a --best --strata -S -p 2 --un ' 
        bowtie_1 += '%s%s_nomatch.fastq --max %s%s_multi.fastq --al %s%s_match.fastq %s ' 
        bowtie_1 += '%s%s %s%s 1>>%s 2>%s' 
        
        # bowtie will only add info to log
        bowtie =  '%s -v 2 -y -m 1 -a --best --strata -S -p 2 --un ' 
        bowtie += '%s%s_nomatch.fastq --max %s%s_multi.fastq --al %s%s_match.fastq %s ' 
        bowtie += '%s%s %s%s 1>>%s 2>>%s' 
        
        # first, align to ladder index to subtract
        bowtie_ladder = bowtie_1 %   (paths_in['path_bowtie'], paths_out['path_ladder'], fname, 
                                         paths_out['path_ladder'], fname, paths_out['path_ladder'], fname, 
                                         paths_in['btindex_ladder'], 
                                         paths_out['path_filter'], fname + '-trimmed.fastq', 
                                         paths_out['path_temp'], fname + '_ladder_match.SAM',
                                         file_log, file_log)
        ladder.append(bowtie_ladder)

        # second, align to ladder index to subtract
        bowtie_tRNA = bowtie %     (paths_in['path_bowtie'], paths_out['path_trna'], fname, 
                                         paths_out['path_trna'], fname, paths_out['path_trna'], fname, 
                                         paths_in['btindex_trna'], 
                                         paths_out['path_ladder'], fname + '_nomatch.fastq', 
                                         paths_out['path_temp'], fname + '_tRNA_match.SAM',
                                         file_log, file_log)
        tRNA.append(bowtie_tRNA)

        # third, align to the rRNA index        
        bowtie_rRNA = bowtie %     (paths_in['path_bowtie'], paths_out['path_rrna'], fname, 
                                         paths_out['path_rrna'], fname, paths_out['path_rrna'], fname, 
                                         paths_in['btindex_rrna'], 
                                         paths_out['path_trna'], fname + '_nomatch.fastq', 
                                         paths_out['path_temp'], fname + '_rRNA_match.SAM',
                                         file_log, file_log)
        rRNA.append(bowtie_rRNA)

        # then align to the chr index  
        bowtie_chr = bowtie %      (paths_in['path_bowtie'], paths_out['path_chr'], fname, 
                                         paths_out['path_chr'], fname, paths_out['path_chr'], fname, 
                                         paths_in['btindex_chr'], 
                                         paths_out['path_rrna'], fname + '_nomatch.fastq', 
                                         paths_out['path_chr'], fname + '_match.SAM',
                                         file_log, file_log)
        chromosome.append(bowtie_chr)

    print "\n------ALIGN------"
    print '\nFiles to align: ' + ', '.join(files)
    print "\n\tStarted Bowtie alignment at " + str(datetime.now())
            
    ribo_util.subprocess_wf(ladder, threads)
    print "\tFinished ladder removal at " + str(datetime.now())

    ribo_util.subprocess_wf(tRNA, threads)
    print "\tFinished tRNA removal at " + str(datetime.now())

    ribo_util.subprocess_wf(rRNA, threads)
    print "\tFinished rRNA removal at " + str(datetime.now())

    ribo_util.subprocess_wf(chromosome, threads)
    print "\tFinished chromosome alignment at " + str(datetime.now())
    
    print "\tCOMPLETED ALIGNING"
    
    return
        
        
        

############################
#####      Density     #####   
############################


def density_3(fname, chr_sam, minlength, maxlength, path_wig, path_den, path_gff, data_type):
    '''Density will be a size separated dictionary = {length : [reads at 0, reads at 1, ....]}
        this makes it easier to select a size range later for analysis'''    
    
    fname = fname
    chr_sam = chr_sam
    minlength = minlength
    maxlength = maxlength
    GFFgen = GFF.parse(path_gff)
    data_type = data_type
    
    if data_type == 'ribo_seq':
        plus_strand  = '0'
        minus_strand = '16'
    else:
        plus_strand  = '16'
        minus_strand = '0'

    # open chr aligned sam file
    f_samfile = open(chr_sam)
    samfile = csv.reader(f_samfile,delimiter='	')
    
    # dictionaries to hold read counts
    density_plus = {}
    density_minus = {}
    density_plus_sizesep = {}
    density_minus_sizesep = {}
    
    if minlength < 0 or maxlength < 0:
        print "Error. Length input not valid."
        return(0)
    
    # Makes 2 sets of indices, one for all reads, and another for size separated:
    for sequence in GFFgen:
        density_plus[sequence.id]  = [0 for x in range(len(sequence))]
        density_minus[sequence.id] = [0 for x in range(len(sequence))]
    
    for length in range(minlength, maxlength + 1):
        density_plus_sizesep[length]  = [0 for x in range(len(sequence))]
        density_minus_sizesep[length] = [0 for x in range(len(sequence))]
        
    total_reads = 0
    mapped_reads = 0

    # Loop through the samfile.
    for read in samfile:
        if read[0][0] == '@':   # Ignore header lines.
            continue

        if read[1] == '4':      # A bowtie mismatch.  
            continue

        chrom = read[2]             # chromosome identified for read in bowtie
        readid = read[0]            # read id
        startp = int(read[3]) -1    # start position. Need to subtract 1 since genomic sequence starts at 1,
        seq = Seq.Seq(read[9])      # sequence of the read
        length = len(seq)           # length of read
        
        if chrom not in density_plus.keys():
            print "Error: Bowtie index and GFF do not match"
        
        total_reads += 1

        # Note that Bowtie reverse complements any sequence aligning to the reverse strand.  
        # and so read[3] is the 3'-end of minus strand reads 

        # Filter to get rid of reads of particular length
        if (length < minlength or length > maxlength):
            continue

        mapped_reads += 1

        # 16 is the minus strand, 0 is the plus strand
        if (read[1] == minus_strand):
            start = startp
            density_minus[chrom][start] += 1 
            density_minus_sizesep[length][start] += 1 

        if (read[1] == plus_strand):
            start = startp + length - 1
            density_plus[chrom][start] += 1 
            density_plus_sizesep[length][start] += 1 
    
    path_oldformat = path_den+"binary/"
    if not os.path.exists(path_oldformat):
            os.makedirs(path_oldformat)
    
    density_plus[sequence.id]  = [float(i) * 1000000 / float(mapped_reads) for i in density_plus[sequence.id]]
    density_minus[sequence.id] = [float(i) * 1000000 / float(mapped_reads) for i in density_minus[sequence.id]]
    
    ribo_util.writebin(density_plus,path_oldformat+fname+"_plus_")
    ribo_util.makePickle(density_plus,path_den+"plus")
    ribo_util.makePickle(density_plus_sizesep,path_den+"plus_sizesep")
    ribo_util.countstowig(density_plus,path_wig+"_plus")
    
    ribo_util.writebin(density_minus,path_oldformat+fname+"_minus_")
    ribo_util.makePickle(density_minus,path_den+"minus")
    ribo_util.makePickle(density_minus_sizesep,path_den+"minus_sizesep")
    ribo_util.countstowig(density_minus,path_wig+"_minus")        


def run_density(inputs, paths_in, paths_out): # all arguments = dict
    
    files     = inputs['files']
    run       = inputs['run_density']
    minlength = inputs['minlength']
    maxlength = inputs['maxlength']
    threads   = inputs['threads'] 
    data_type = inputs['data_type']
    
    if not files:
        print("There are no files")
        return
    
    print "\n-----DENSITY-----"
    print '\nFiles to condense: ' + ', '.join(files)
    print "\n\tStarted density at " + str(datetime.now())
    
    arguments = []
    
    for fname in files:

        # make paths for density files
          
        path_d = paths_out['path_density'] + fname + "/" 
        path_w = paths_out['path_density'] + "wigfiles/" 
        if not os.path.exists(path_d):
            os.makedirs(path_d)
        if not os.path.exists(path_w):
            os.makedirs(path_w)	

        path_den = path_d
        path_wig = path_w + fname
        path_sam = paths_out['path_chr'] + fname + '_match.SAM'
        path_gff = paths_in['path_gff']
        
        if not run == 'yes':
            if not os.path.exists(path_den+"plus"):
                print "ERROR: " + fname + " has no density, change run settings"
                continue
            else: 
                print fname + " has density file"
                continue
                
        if not os.path.exists(path_sam):
            print "ERROR: " + fname + " has no alignment file, has been removed from analysis"
            inputs['files'].remove(fname)
            continue
       
        argument = [fname, path_sam, minlength, maxlength, path_wig, path_den, path_gff, data_type]
        arguments.append(argument)
    
    ribo_util.multiprocess(density_3, arguments, threads)
            
    
    print "\tFinished density at " + str(datetime.now())

    print "\tCOMPLETED DENSITY"
    
##########################
    
def density_adjusted(fname, chr_sam, minlength, maxlength, path_wig, path_den, path_gff):
    '''Density will be a size separated dictionary = {length : [reads at 0, reads at 1, ....]}
        this makes it easier to select a size range later for analysis
        
        adjusted: will shift reads larger than 24 to alignn 3' end'''    
    
    fname = fname
    chr_sam = chr_sam
    minlength = minlength
    maxlength = maxlength
    GFFgen = GFF.parse(path_gff)

    # open chr aligned sam file
    f_samfile = open(chr_sam)
    samfile = csv.reader(f_samfile,delimiter='	')
    
    # dictionaries to hold read counts
    density_plus = {}
    density_minus = {}
    density_plus_sizesep = {}
    density_minus_sizesep = {}
    
    if minlength < 0 or maxlength < 0:
        print "Error. Length input not valid."
        return(0)
    
    # Makes 2 sets of indices, one for all reads, and another for size separated:
    for sequence in GFFgen:
        density_plus[sequence.id]  = [0 for x in range(len(sequence)+20)]
        density_minus[sequence.id] = [0 for x in range(len(sequence)+20)]
    
    for length in range(minlength, maxlength + 1):
        density_plus_sizesep[length]  = [0 for x in range(len(sequence)+20)]
        density_minus_sizesep[length] = [0 for x in range(len(sequence)+20)]

    total_reads = 0
    mapped_reads = 0

    # Loop through the samfile.
    for read in samfile:
        if read[0][0] == '@':   # Ignore header lines.
            continue

        if read[1] == '4':      # A bowtie mismatch.  
            continue

        chrom = read[2]             # chromosome identified for read in bowtie
        readid = read[0]            # read id
        startp = int(read[3]) -1    # start position. Need to subtract 1 since genomic sequence starts at 1,
        seq = Seq.Seq(read[9])      # sequence of the read
        length = len(seq)           # length of read
        
        if length < 23:
            length_shift = 24 - length
        else:
            length_shift = 0
        
        if chrom not in density_plus.keys():
            print "Error: Bowtie index and GFF do not match"
        
        total_reads += 1

        # Note that Bowtie reverse complements any sequence aligning to the reverse strand.  
        # and so read[3] is the 3'-end of minus strand reads 

        # Filter to get rid of reads of particular length. Or a particular strand.
        if (length < minlength or length > maxlength):
            continue

        mapped_reads += 1

        # 16 is the minus strand, 0 is the plus strand
        if (read[1] == '16'):
            start = startp - length_shift
            density_minus[chrom][start] += 1 
            density_minus_sizesep[length][start] += 1 

        if (read[1] == '0'):
            start = startp + length - 1 + length_shift
            density_plus[chrom][start] += 1 
            density_plus_sizesep[length][start] += 1 
    
    path_oldformat = path_den+"binary/"
    if not os.path.exists(path_oldformat):
            os.makedirs(path_oldformat)
    
    density_plus[sequence.id]  = [float(i) * 1000000 / float(mapped_reads) for i in density_plus[sequence.id]]
    density_minus[sequence.id] = [float(i) * 1000000 / float(mapped_reads) for i in density_minus[sequence.id]]
    
    ribo_util.writebin(density_plus,path_oldformat+fname+"_plus_")
    ribo_util.makePickle(density_plus,path_den+"plus")
    ribo_util.makePickle(density_plus_sizesep,path_den+"plus_sizesep")
    ribo_util.countstowig(density_plus,path_wig+"_plus")
    
    ribo_util.writebin(density_minus,path_oldformat+fname+"_minus_")
    ribo_util.makePickle(density_minus,path_den+"minus")
    ribo_util.makePickle(density_minus_sizesep,path_den+"minus_sizesep")
    ribo_util.countstowig(density_minus,path_wig+"_minus")        


def run_density_adjusted(inputs, paths_in, paths_out): # all arguments = dict
    
    files     = inputs['files']
    run       = inputs['run_density']
    minlength = inputs['minlength']
    maxlength = inputs['maxlength']
    threads   = inputs['threads'] 
    
    if not files:
        print("There are no files")
        return
    
    print "\n-----DENSITY-----"
    print '\nFiles to condense: ' + ', '.join(files)
    print "\n\tStarted density at " + str(datetime.now())
    
    arguments = []
    
    for fname in files:

        # make paths for density files
          
        path_d = paths_out['path_density'] + fname + "/adjusted/" 
        path_w = paths_out['path_density'] + "wigfiles/adjusted/" 
        if not os.path.exists(path_d):
            os.makedirs(path_d)
        if not os.path.exists(path_w):
            os.makedirs(path_w)	

        path_den = path_d
        path_wig = path_w + fname
        path_sam = paths_out['path_chr'] + fname + '_match.SAM'
        path_gff = paths_in['path_gff']
        
        if not run == 'yes':
            if not os.path.exists(path_den+"plus"):
                print "ERROR: " + fname + " has no density, change run settings"
                continue
            else: 
                print fname + " has density file"
                continue
                
        if not os.path.exists(path_sam):
            print "ERROR: " + fname + " has no alignment file, has been removed from analysis"
            inputs['files'].remove(fname)
            continue
       
        argument = [fname, path_sam, minlength, maxlength, path_wig, path_den, path_gff]
        arguments.append(argument)
    
    ribo_util.multiprocess(density_adjusted, arguments, threads)
            
    
    print "\tFinished density at " + str(datetime.now())

    print "\tCOMPLETED DENSITY"