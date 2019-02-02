'''
    Misc utilities to process and analyze profiling data
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
from multiprocessing import Process, Pool
import os, time
import subprocess
import struct
import cPickle as pickle
import csv
from BCBio import GFF
from Bio import Seq
import itertools
import sys

import pandas as pd
import numpy as np
from IPython.display import display
from collections import deque


'''
Table of Contents:

    -- SYSTEM UTILITIES:
    
        - multiprocess      = running multiple samples through the same function in parallel
        - subprocess_wf     = running multiple os commands in parallel
        - createpath        = make file in indicated path
        - check_inputs      = query if files exist
        - rename_FASTQ      = rename and concat Fastq files from sequencing
        
    -- PROFILING UTILITIES:
    
        - log
        - get_log
        - GFF_to_dict
        - nextgene
        - merge_density_lenghts
        - getallcounts
        - get_density_rpm
        - getRPKM
        - get_genetic_code
        - orf_motif_positions
        - 3to5_alignment
        
    -- DATA STROAGE AND MANIPULATION:
    
        - writebin
        - countstowig
        - makePickle
        - unPickle
        - loadlargePickles
        - GFF_to_dict
        - heatmapdict_to_df
        - dict_to_df
    
        
'''

############################
##### SYSTEM UTILITIES #####
############################



def multiprocess(function, arguments, threads):
    ''' 
    Multiprocess utility to run multiple functions in parallel
    Useful for creating density or running analysis on multiple files
    '''
    processes = []
    
    for arg in arguments:
        p = Process(target = function, args = arg)
        processes.append(p)
        p.start()
    
    for p in processes:
        p.join()
        
    return
    

def subprocess_wf(cmds, threads):
    ''' 
    Subprocess utility to run subprocess commands in parallel
    Useful for filtering or bowtie when runnning multiple samples at once
    Skewer uses multiple threads, best to use 2 threads if hyperthreading, else 1
    Bowtie uses 2 threads, best to use 4 for hyperthreading, else 2
    '''
    
    if not cmds: return # empty list
 
    def done(p):
        return p.poll() is not None
    def success(p):
        return p.returncode == 0
    def fail():
        sys.exit(1)

    threads = threads
    
    processes = []
    while True:
        while cmds and len(processes) <= threads:
            task = cmds.pop()
            processes.append(subprocess.Popen(task, shell=True))

        for p in processes:
            if done(p):
                if success(p):
                    processes.remove(p)
                else:
                    fail()

        if not processes and not cmds:
            break
        else:
            time.sleep(0.05)

            
def createpath(inputs, paths_dict):
    '''
    Creates paths for a dictionary of path values if they dont exist
    '''
    for path in paths_dict.values(): 
        if not os.path.exists(path):
            os.makedirs(path)
    
    '''
    For path_analysis, make folder with fname
    '''
    files = inputs['files']

    if 'path_analysis' in paths_dict:
        for fname in files:
            path_analysis = paths_dict['path_analysis'] + fname + "/" 
            
            if not os.path.exists(path_analysis):
                os.makedirs(path_analysis)
                
    if 'path_figures' in paths_dict:
        for fname in files:
            path_figures = paths_dict['path_figures'] + fname + "/" 
            
            if not os.path.exists(path_figures):
                os.makedirs(path_figures)
        
    
    ### for ribosome profiling ###
    # Output structure =  ...data   /libraries/all_libraries/FASTQ              /fname
    #                               /gff                    /filterdata         /fname
    #                               /scripts                /alignments
    #                                                           /chr            /fname
    #                                                           /tRNA           /fname
    #                                                           /rRNA           /fname
    #                                                       /density
    #                                                       /all            /fname
    #                                                       /separated      /fname_length
    #                                                       /tmpds
    #

    
def check_inputs(inputs, paths_in, step):
    
    files = inputs['files']
    no_file = []
    
    if step == 'density':
        for fname in files: 
            if not os.path.exists(paths_in['path_fastq'] + fname):
                no_file.append(fname)
    else: 
        for fname in files: 
            if not os.path.exists(paths_in['path_density'] + fname):
                no_file.append(fname)
    
    for fname in no_file:
        files.remove(fname)
        print fname + ' does not have a FASTQ/density file and was removed'
    
    for path_name, path in paths_in.iteritems():
        if 'indexes' in path:
            path = path + '.1.ebwt'
        if not os.path.exists(path):
            print 'The path for ' + path_name + ' does not exist'
            

def rename_FASTQ(inputs, library_creator, paths_in, paths_out):
    
    '''Function will take data downloaded from GRCF (placed in order_path), concat it if 
    multiple files exist, and move/rename it
    
    Make sure order form uses FM## naming (SM_Tag column on order sheet)!!
    '''
    
    files          = inputs['files']
    order_name     = inputs['order_name']
    order_path     = paths_in['fastq_download'] + order_name + '/'
    order_csv      = order_path + order_name + '.csv'
    order_sheet_df = pd.read_csv(order_csv)
    
    rename_path  = paths_in['path_fastq']
    
    '''Order Sheet Layout: 

            Project	            FCID	  Lane	Index	SM_Tag	File_Name
        abuskir1_FM_EL_143379	HCJWKBCX2	1	GCCAAT	FM80	HCJWKBCX2_1_GCCAAT
        abuskir1_FM_EL_143379	HCJWKBCX2	1	GTCCGC	FM83	HCJWKBCX2_1_GTCCGC
        abuskir1_FM_EL_143379	HCJWKBCX2	1	ACAGTG	FM79	HCJWKBCX2_1_ACAGTG'''
        
    # skip renaming if not needed/download file doesnt exist
    if order_name == 'none':
        print 'Skipped FASTQ rename'
        return
    elif not os.path.exists(order_path):
        print 'Download folder does not exist. Skipped FASTQ rename'
        return
    
    print "\nStarted Fastq renaming at " + str(datetime.now())
    
    # iterate through files:
    for fname in files:
        
        rename_file = rename_path + fname
        
        # ignore fname if already renamed
        if os.path.exists(rename_file):
            print fname + ' has been renamed'
            continue
            
        # find all rows in order sheet with fname        
        fname_df = order_sheet_df.query('SM_Tag == "' + fname + '"') 
        
        #if there are multiple files to concatonate
        if fname_df.shape[0] > 1: 
            
            # iterate through order_sheet to find file info
            sample_paths = ''
            for index, row in fname_df.iterrows():
                sample_name = row['File_Name'] + '_1.fastq.gz'
                sample_path = order_path + sample_name
                                
                if not os.path.exists(sample_path):
                    print sample_name + ' missing'
                    continue
                else: 
                    sample_paths += sample_path + ' '
            
            print fname + ' has multuple files, concatonating and renaming' 
           
            # combine files using os command
            concat_gz_out  = order_path + fname + '.gz'     
            concat_command = 'cat %s> %s' % (sample_paths, concat_gz_out)
            
            os.system(concat_command)
            
            # unzip gz file
            os.system('gunzip ' + concat_gz_out)
            unzip_out = order_path + fname
            
            # move file to destination
            os.system('mv ' + unzip_out + ' ' + rename_path)
        
        # if there is only one file:
        else: 
            for index, row in fname_df.iterrows():
                sample_name = row['File_Name'] + '_1.fastq.gz'
                sample_path = order_path + sample_name
                
                
                if not os.path.exists(sample_path):
                    print sample_name + ' missing'
                    continue
                
                print fname + ' has a single file, renaming' 

                # uzip gz file:
                os.system('gunzip ' + sample_path)
                unzip_out = order_path + row['File_Name'] + '_1.fastq'
                
                 # move and rename file to destination
                os.system('mv ' + unzip_out + ' ' + rename_file)
    
    print "\n\tFinished Fastq renaming at " + str(datetime.now())

    
                
            
###############################            
##### PROFILING UTILITIES #####
###############################



def log(fname, paths_out, text):
    '''
    Adds text to log
    '''
    
    f = open(paths_out['path_log'] + fname, 'a')
    f.write(text) 
    f.close()
    
    
def get_filter_bowtie_log(inputs, paths_in, paths_out):
    '''
    Reads log files to extract filter and bowtie stats
    '''
    
    files = inputs['files']
    
    if not files:
        print("There are no files")
        return
    
    # processes are the row indeces for pandas
    processes = [      
        'total_reads',
        'filtered_reads',
        'ladder',
        'tRNA',
        'rRNA',
        'chromosome',
        '% Total mapped',
        'unaligned'
        ]
    
    #processing_data contains log info for all files
    processing_data = {}
    log_data        = {}
    for fname in files:
        
        processing_data[fname] = []
        bowtie_data = [0] * 12   
        
        for process in ['_filter', '_bowtie']:
            f = open(paths_out['path_log'] + fname + process)
        
            #bowtie_data contains log info for each file
            counter = 0 
            
            for line in f:
                if line.strip().endswith('reads processed; of these:'):
                    total_reads = int(line.strip().split()[0])
                    processing_data[fname].append(total_reads)
                if line.strip().endswith('reads available; of these:'):
                    filtered_reads = int(line.strip().split()[0])     # reads that passed filtering
                    processing_data[fname].append(filtered_reads)
                if line.strip().startswith('# reads with at'): 
                    aligned = int(line.strip().split()[-2])
                    if counter in [0, 3, 6, 9]: 
                        bowtie_data[counter] += aligned
                        counter += 1
                    else:
                        counter += 1
                        bowtie_data[counter] += aligned
                if line.strip().startswith('# reads that failed'):
                    unaligned =  int(line.strip().split()[-2])
                    if counter in [1, 4, 7, 10]: 
                        bowtie_data[counter] += unaligned
                        counter += 1
                    else:
                        counter += 1
                        bowtie_data[counter] += unaligned
                if line.strip().startswith('# reads with alignments'):
                    multiple = int(line.strip().split()[-2])
                    if counter in [2, 5, 8, 11]: 
                        bowtie_data[counter] += multiple
                        counter += 1
                    else:
                        counter += 1
                        bowtie_data[counter] += multiple
                if line.strip().startswith('No alignments'):
                    No_alignments = 0
                    bowtie_data.append(No_alignments) 
        
        # do the math
        reads_filtered = total_reads - filtered_reads   # reads removed by filtering
        ladder_reads   = bowtie_data[0] + bowtie_data[2]
        tRNA_reads     = bowtie_data[3] + bowtie_data[5]
        rRNA_reads     = bowtie_data[6] + bowtie_data[8]
        chr_single     = bowtie_data[9]
        percent_mapped = int((float(bowtie_data[9]) / total_reads) * 100)
        no_alignment   = bowtie_data[10] + bowtie_data[11]
        
        processing_data[fname].append(ladder_reads)
        processing_data[fname].append(tRNA_reads)
        processing_data[fname].append(rRNA_reads)
        processing_data[fname].append(chr_single)
        processing_data[fname].append(percent_mapped)
        processing_data[fname].append(no_alignment)
        
        f.close()
        
        #Add skewer and bowtie output to analysis log:
        
        log_data['analysis_breakdown'] = {'Total Reads': total_reads, 'Reads Filtered': reads_filtered,
                                            'Ladder': ladder_reads, 'tRNA': tRNA_reads, 'rRNA': rRNA_reads, 
                                            'Chromosome': chr_single, 'Other': no_alignment}
        
        log_function = 'ribo_density'
        analysis_log(fname, log_function, log_data, paths_in, paths_out)
        
        # make dataframe to view
        df = pd.DataFrame(data = processing_data ,index=processes)
    return display(df)

def analysis_log(fname, log_function, log_data, paths_in, paths_out):
    '''
    Stores settings information for samples
    
    log_function = name of function being executed, eg: ribo_avggene
    log_data     = {settings : []} or {analysis_breakdown : []}
    
    log_dictionary = {  ribo_density : { settings           : { linker_seq, cutoff, minlength, maxlength} 
                                         analysis_breakdown : {total_reads, filtered reads, mapped_ladder, 
                                                              mapped_tRNA, mapped_rRNA, mapped_chr, unaligned} 
                                                             
                        ribo_avggene : { settings    : {'size range', 'length in orf', 'length out orf', 
                                                           'density type', 'weighting'], 
                                        analysis_breakdown :[genes_in_avg]
                                                            nextgene_cut : int
                                                            length_cut   : int
                                                            RPM_cut      : int}, 
                                        
                        etc... 
    '''
    log_path = paths_out['path_analysis_log'] + fname + '/'
    log_file = log_path + fname
    
    if not os.path.exists(log_file):
        os.makedirs(log_path)
        log_all_data = {}
    else:
        log_all_data = unPickle(log_file)
    
    # add log_data to log file:
    log_all_data[log_function] = log_data
    makePickle(log_all_data, log_file)
    

def GFF_to_dict(paths_in, gff_settings):
    
    '''Parse gff into dict:
        - feat_of_interest = what to look for in gff (protein_coding, tRNA, rRNA, etc)
        - name_qual        = qualifier for alias/gene name (Name, gene_id)
        - name_qual_alt    = alternative qualifier, if none, set as 'none' 
        - biotype_qual     = qualifier for type of feature (biotype, etc)
        
        These values must correspont to values in the GFF'''
    
    '''Unload gff_settings'''
    
    path_out         = gff_settings['path_out']
    feat_of_interest = gff_settings['feat_of_interest']  #all, protein_coding, tRNA, rRNA
    name_qual        = gff_settings['name_qual']
    name_qual_alt    = gff_settings['name_qual_alt']
    remove_genes     = gff_settings['remove_genes']
    aSD_seq          = gff_settings['aSD_seq']
    path_badgenes    = paths_in['path_badgenes']
    
    '''Output path can be defined, or use 0 to set as the annotation file for my main pipeline'''
    
    if path_out == 0:
        path_gff_dict = paths_in['path_gff_dict']
    else:
        path_gff_dict = path_out
    
    '''Parse GFF using BCBio'''
    
    GFFgen = GFF.parse(paths_in['path_gff'])
    chr = GFFgen.next()
    feat_num = 0
    
    '''Define genetic code'''
    
    aa_code, codon_code = ribo_util.get_genetic_code()
    aa_comp_dict = {}
    
    '''Define data arrays: will be used as columns for pandas DateFrame'''

    gff_dict   = {}
    aliaslist  = []
    startlist  = []
    stoplist   = []
    seqlist    = []
    typelist   = []
    strandlist = []
    startcodon = []
    stopcodon  = []
    SDaffinity = []
    IEpoint    = []
    G_content  = []
    C_content  = []
    A_content  = []
    T_content  = []
    
    aa_start = {}
    aa_stop  = {}
    
    for aa in aa_code.keys():
        aa_start[aa] = [0] * 25
        aa_stop[aa]  = [0] * 25
    
    '''Make list of bad genes'''
            

    # from Gene-Wei-Li 

    bad_genes = pd.read_csv(path_badgenes)
    bad_genes = bad_genes.to_dict(orient='list')
    bad_genes = bad_genes['GeneName']
        
        

    '''Sift through GFF for relevant information'''
    
    for feature in chr.features:
        
        if feature.sub_features == []:
                feat_num+=1
                continue
                
        if remove_genes == 'yes':
            
            '''Skip over non-CDS annotations'''
            
            if not feature.sub_features[0].type == feat_of_interest:
                feat_num+=1
                continue        
            elif feature.qualifiers.has_key('pseudo') == True:
                feat_num+=1
                continue
            else:
                feature_type = 'CDS'
        else: 
            
            '''Add feat type to GFF, noting pseudogenes'''
            
            if feature.qualifiers.has_key('pseudo') == True:
                feature_type = 'pseudo'
            else: 
                feature_type = feature.sub_features[0].type 
                
                    
        '''Get feature name'''
        
        if name_qual in feature.qualifiers:
            feat_name = feature.qualifiers[name_qual][0]
        elif name_qual_alt in feature.qualifiers:
            feat_name = feature.qualifiers[name_qual_alt][0]
        else:
            feat_name = 'None'
            feat_num+=1
            continue
        
        '''Remove feature if bad'''
        
        if remove_genes == 'yes':
            if feat_name in bad_genes:
                feat_num+=1
                continue
        else: 
            if feat_name in bad_genes:
                feature_type = 'bad'
                            
        
        '''Get start, end, and strand position'''
        
        start  = feature.location.start.position 
        end    = feature.location.end.position
        strand = feature.strand           
       
        '''Analyze features of interest for feat information'''   
        

        alias = feat_name
        
            
        '''Each strand is treated differently, + strand == 1'''
            
        if strand == 1:
                
            '''I save gene sequence + 50 bp from each end:
            makes it easier to analyze start and stop sequence 
            context without using whole genome sequence'''

            if start < 50:      # if gene is near the beginning of genome sequence:
                sequence = 'N' * (50 - start)                  # TB GFF starts at 0, add N * 50
                sequence = sequence + chr[0:end+50].seq        # gene sequence + 50nt at each end
            else: 
                sequence = chr[start-50:end+50].seq            # gene sequence + 50nt at each end

            strand_val = '+'
            startcodon_pos = start
            stopcodon_pos  = end-1

            if start > 200: 
                upstream_seq = chr[start-200:start+100].seq

        else:

            '''For minus strand, 'end' is start codon, 'start' is stop codon
            and sequence is reverse compliment of gene sequence.'''

            sequence_rc = chr[start-50:end+50].seq
            sequence = sequence_rc.reverse_complement()

            strand_val = '-'
            startcodon_pos = end-1
            stopcodon_pos  = start

            if end + 200 > len(chr.seq):
                upstream_seq = 'none'
            else:
                upstream_seq_rc = chr[end-100:end+200].seq
                upstream_seq = upstream_seq_rc.reverse_complement()
            
        sequence = str(sequence)
        start_codon = sequence[50:53:1]
        stop_codon  = sequence[-53:-50]
        
        '''print an alias and sequence of interest'''
        
        if alias == 'raiA':
            
            print alias, sequence 

        '''get sequence from start to stop for GC analysis'''

        CDS_seq = sequence[50:-50:1]
        CDS_seq_length = len(CDS_seq)

        G, C, A, T = GC_of_CDS(CDS_seq)
        
        '''get codon composition metaplot data'''
                
        if CDS_seq_length % 3 == 0: 
            codons_seq = [CDS_seq[i:i+3] for i in range(0, CDS_seq_length, 3)]
            aa_seq     = [codon_code[codon] for codon in codons_seq] 
        else: 
            aa_seq = []

        #get start and stop codon alignment
        
        if len(aa_seq) > 50:
            
            position_index = range(0, 25)
            for position in position_index:

                # start position data
                aa = aa_seq[position]
                aa_start[aa][position] += 1

                # stop position data
                aa = aa_seq[-25 + position]
                aa_stop[aa][position] += 1  

        '''Calculate SD affinity'''

        SD_seq = sequence[30:50:1]    # analyze 20 nt upstream of start codons
        SD_affinity = shine_dalgarno_affinity(aSD_seq, SD_seq)

        '''Calculate Isoelectric point of protein'''
        
        analysed_seq = ProteinAnalysis(CDS_seq)
        pI           = analysed_seq.isoelectric_point()
        
        '''Append data to lists'''
            
        typelist.append(feature_type)
        aliaslist.append(alias)
        seqlist.append(sequence)
        strandlist.append(strand_val)
        startlist.append(startcodon_pos)
        stoplist.append(stopcodon_pos)
        startcodon.append(start_codon)
        stopcodon.append(stop_codon)
        SDaffinity.append(SD_affinity)
        IEpoint.append(pI)
        G_content.append(G)
        C_content.append(C)
        A_content.append(A)
        T_content.append(T)

        feat_num+=1

    '''Normalize codon composition data'''
    
    for aa in aa_start.keys():
        total_aa_start = sum(aa_start[aa])
        avg_aa_start   = float(total_aa_start / 25)
        if avg_aa_start == 0:
            avg_aa_start = 1
        aa_start[aa]   = [ float(i / avg_aa_start) for i in aa_start[aa]]
        
        total_aa_stop = sum(aa_stop[aa])
        avg_aa_stop   = float(total_aa_stop / 25)
        if avg_aa_stop == 0:
            avg_aa_stop = 1
        aa_stop[aa]   = [ float(i / avg_aa_stop) for i in aa_stop[aa]]
    
    # reorder data for plotting
    data_df = {}
    position_list    = []
    aa_list          = []
    value_list_start = []
    value_list_stop  = []
    aa_property      = []
    
    
    
    for position in range(0, 25):
        for aa in aa_start.keys():
            
            if aa in ['R', 'H', 'K']:
                aa_prop = 'positive'
            if aa in ['D', 'E']:
                aa_prop = 'negative'
            if aa in ['S', 'T', 'N', 'Q']:
                aa_prop = 'polar'
            if aa in ['A', 'V', 'I', 'L', 'M', 'C']:
                aa_prop = 'hydrophobic'
            if aa in ['G', 'P']:
                aa_prop = 'steric'
            if aa in ['F', 'Y', 'W']:
                aa_prop = 'bulky'
            if aa in ['_']:
                continue
            
            value_start = aa_start[aa][position]
            value_stop  = aa_stop[aa][position]
            value_list_start.append(value_start)
            value_list_stop.append(value_stop)
            position_list.append(position)
            aa_list.append(aa)
            aa_property.append(aa_prop)
            
    data_df['Start']      = value_list_start
    data_df['Stop']       = value_list_stop
    data_df['Amino_Acid'] = aa_list
    data_df['Position']   = position_list
    data_df['Property']   = aa_property
    
    data_df = pd.DataFrame(data_df)
    
    display(data_df)
    
    sns.set(style="white")    
    plt.figure(figsize=(25,6))
    plt.subplot(1,2,1)
    plot1 = sns.lineplot(x="Position", y="Start", hue="Amino_Acid", data=data_df)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.subplot(1,2,2)
    plot2 = sns.lineplot(x="Position", y="Stop", hue="Amino_Acid", data=data_df)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.show()
    
    
    '''Append lists to gff_dict'''

    gff_dict['Alias']    = aliaslist
    gff_dict['Strand']   = strandlist
    gff_dict['Start']    = startlist
    gff_dict['Stop']     = stoplist
    gff_dict['Sequence'] = seqlist
    gff_dict['Start_Codon'] = startcodon
    gff_dict['Stop_Codon']  = stopcodon
    gff_dict['Type']        = typelist
    gff_dict['SD_affinity'] = SDaffinity
    gff_dict['pI']        = SDaffinity
    gff_dict['G_content'] = G_content
    gff_dict['C_content'] = C_content
    gff_dict['A_content'] = A_content
    gff_dict['T_content'] = T_content

           
    '''Pickle dict for use later'''
    ribo_util.makePickle(gff_dict,path_gff_dict)
    
    '''print dataframe, and save as .csv for use later'''
    ## Print GFF to check 
    gff_df = pd.DataFrame(gff_dict)
    display(gff_df)
    gff_df.to_csv(path_gff_dict + '.csv')
    
    return data_df


def nextgene(alias, gff_dict):
    '''find the next gene in the strand'''
    
    alias_list  = gff_dict['Alias'] 
    strand_list = gff_dict['Strand'] 
    start_list  = gff_dict['Start'] 
    stop_list   = gff_dict['Stop']
    
    index  = alias_list.index(alias)
    strand = strand_list[index]
    maxindex = len(alias_list) - 1
    nextgene = {}
    
    if strand == '+':
        next_index = index + 1
        
        if next_index > maxindex:                #if last gene in '+' strand
            nextgene['alias']    = 'none'
            nextgene['strand']   = ''
            nextgene['start']    = ''
            nextgene['distance'] = 10000000
            nextgene['stop']     = ''
            return nextgene
            
        while strand_list[next_index] != '+':    #skips through '-' strand genesif they are next
            next_index += 1 
            if next_index > maxindex:            #if last gene in '+'
                nextgene['alias']    = 'none'
                nextgene['strand']   = ''
                nextgene['start']    = ''
                nextgene['distance'] = 10000000
                nextgene['stop']     = ''
                return nextgene
        
        distance = start_list[next_index] - stop_list[index] # get distance between genes
        
    if strand == '-':
        next_index = index - 1
        
        if next_index < 0:
            nextgene['alias']    = 'none'
            nextgene['strand']   = ''
            nextgene['start']    = ''
            nextgene['distance'] = 10000000
            nextgene['stop']     = ''
            return nextgene
            
        while strand_list[next_index] != '-':
            next_index -= 1 
            if next_index < 0:
                nextgene['alias']    = 'none'
                nextgene['strand']   = ''
                nextgene['start']    = ''
                nextgene['distance'] = 10000000
                nextgene['stop']     = ''
                return nextgene
            
        distance = stop_list[index] - start_list[next_index] 
    
    nextgene['alias']    = alias_list[next_index]
    nextgene['strand']   = strand_list[next_index]
    nextgene['distance'] = distance
    nextgene['start']    = start_list[next_index]
    nextgene['stop']     = stop_list[next_index]
    
    return nextgene

def prevgene(alias, gff_dict):
    
    alias_list  = gff_dict['Alias'] 
    strand_list = gff_dict['Strand'] 
    start_list  = gff_dict['Start'] 
    stop_list   = gff_dict['Stop']
    
    index  = alias_list.index(alias)
    strand = strand_list[index]
    maxindex = len(alias_list) - 1
    prevgene = {}
    
    if strand == '+':
        prev_index = index - 1
        
        if prev_index < 0:
            prevgene['alias']    = 'none'
            prevgene['strand']   = ''
            prevgene['start']    = ''
            prevgene['distance'] = 10000000
            prevgene['stop']     = ''
            return prevgene
            
        while strand_list[prev_index] != '+':
            prev_index -= 1 
            if prev_index < 0:
                prevgene['alias']    = 'none'
                prevgene['strand']   = ''
                prevgene['start']    = ''
                prevgene['distance'] = 10000000
                prevgene['stop']     = ''
                return prevgene
            
        distance = start_list[index] - stop_list[prev_index] 
        
    if strand == '-':
        prev_index = index + 1
        
        if prev_index > maxindex:                #if last gene in '-' strand
            prevgene['alias']    = 'none'
            prevgene['strand']   = ''
            prevgene['start']    = ''
            prevgene['distance'] = 10000000
            prevgene['stop']     = ''
            return prevgene
            
        while strand_list[prev_index] != '-':    #skips through '+' strand genesif they are next
            prev_index += 1 
            if prev_index > maxindex:            #if last gene in '+'
                prevgene['alias']    = 'none'
                prevgene['strand']   = ''
                prevgene['start']    = ''
                prevgene['distance'] = 10000000
                prevgene['stop']     = ''
                return prevgene
        
        distance = stop_list[prev_index] - start_list[index] # get distance between genes
        
    
    prevgene['alias']    = alias_list[prev_index]
    prevgene['strand']   = strand_list[prev_index]
    prevgene['distance'] = distance
    prevgene['start']    = start_list[prev_index]
    prevgene['stop']     = stop_list[prev_index]
    
    return prevgene

def merge_density_lenghts(density_plus, density_minus, minlength, maxlength):
    
    lengthindex  = range(minlength, maxlength+1)
    
    length_density = len(density_plus[minlength])
    
    merged_density_plus  = [0] * length_density
    merged_density_minus = [0] * length_density
    
    length_density = len(density_plus[minlength])
    
    for position in range(0, length_density):
        for lenght in lengthindex:
            merged_density_plus[position]  += density_plus[length][position]
            merged_density_minus[position] += density_minus[length][position]
            
    return merged_density_plus, merged_density_minus

def get_allcounts(density_plus, density_minus, minlength, maxlength):
    '''get total reads in density within the given readlengths'''
    
    lengthindex  = range(minlength, maxlength+1)
    
    allcounts_plus  = 0 
    allcounts_minus = 0 
    totalcounts     = 0 
    
    for length in lengthindex:
        allcounts_plus  += sum(density_plus[length])
        allcounts_minus += sum(density_minus[length])
            
    totalcounts = allcounts_plus + allcounts_minus   
    
    return totalcounts


def get_density_rpm(density_plus, density_minus, minlength, maxlength):
    '''convert raw reads to reads per million for slected readlenghts, 
        to normalize for read depth'''
    
    totalcounts = get_allcounts(density_plus, density_minus, minlength, maxlength)
    totalcounts = float(totalcounts)
            
    density_plus_rpm  = {}
    density_minus_rpm = {}
    
    lengthindex  = range(minlength, maxlength+1)
    
    for length in lengthindex:
        densityplus_rpm  = [float(i) / totalcounts * 1000000 for i in density_plus[length]] 
        densityminus_rpm = [float(i) / totalcounts * 1000000 for i in density_minus[length]] 
        
        density_plus_rpm[length]  = densityplus_rpm
        density_minus_rpm[length] = densityminus_rpm
        
        
    return density_plus_rpm, density_minus_rpm
        
    
def get_genecounts(start, stop, strand, density_plus, density_minus, minlength, maxlength):
    '''count number of reads on a gene, can be raw reads or rpm depending on input density'''
    
    lengthindex = range(minlength, maxlength + 1)
    genelength  = abs(start - stop) + 1

    counts_gene   = [0] * genelength
    counts_length = {}
    
    if strand == '+':
        for length in lengthindex:
            reads_length = density_plus[length][start: stop + 1] 
            counts_length[length] = reads_length
            counts_gene = [x + y for x, y in itertools.izip(counts_gene, reads_length)]
           
    if strand == '-':
        for length in lengthindex:
            reads_length = density_minus[length][start: stop + 1: -1] 
            counts_length[length] = reads_length
            counts_gene = [x + y for x, y in itertools.izip(counts_gene, reads_length)]
    
    return counts_gene, counts_length

            
def get_RPKM(alias, start, stop, strand, density_plus_rpm, density_minus_rpm, minlength, maxlength):
    'input rpm density ---> output RPKM'
    genelength = abs(start - stop) + 1
        
    rpm_gene, rpm_sep = get_genecounts(start, stop, strand, density_plus_rpm, density_minus_rpm, minlength, maxlength)
    
    rpm_gene = sum(rpm_gene)
    
    rpkm = rpm_gene / genelength * 1000
    
    return  rpkm

def get_rpc(alias, start, stop, strand, density_plus_rpm, density_minus_rpm, minlength, maxlength):
    'input raw reads density ---> output reads per codon'
    genelength = abs(start - stop) + 1
        
    reads_gene, reads_sep = get_genecounts(start, stop, strand, density_plus_rpm, density_minus_rpm, minlength, maxlength)
    
    reads_gene = sum(reads_gene)
    
    codons_in_gene = genelength / 3
    
    rpc = reads_gene / codons_in_gene
    
    return  rpc


def get_genetic_code():
    
    aa_keys = {
        'I' : ['ATA', 'ATC', 'ATT'],
        'M' : ['ATG'],
        'T' : ['ACA', 'ACC', 'ACG', 'ACT'],
        'N' : ['AAC', 'AAT'],
        'K' : ['AAA', 'AAG'],
        'S' : ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'],
        'R' : ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'],
        'L' : ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'],
        'P' : ['CCA', 'CCC', 'CCG', 'CCT'],
        'H' : ['CAC', 'CAT'],
        'Q' : ['CAA', 'CAG'],
        'V' : ['GTA', 'GTC', 'GTG', 'GTT'],
        'A' : ['GCA', 'GCC', 'GCG', 'GCT'],
        'D' : ['GAC', 'GAT'],
        'E' : ['GAA', 'GAG'],
        'G' : ['GGA', 'GGC', 'GGG', 'GGT'],
        'F' : ['TTC', 'TTT'],
        'Y' : ['TAC', 'TAT'],
        'C' : ['TGC', 'TGT'],
        'W' : ['TGG'],
        '_' : ['TAA', 'TAG', 'TGA']
        }
        
        
        
    codon_keys = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
        }
    
    return aa_keys, codon_keys


def orf_motif_position(motif, gff_dict):
    '''uses gff information to obtain positional information of motifs of interest'''
    
    
    
    return 

def alignment_3to5(density_plus, density_minus, minlength, maxlength):
    '''convert 3' density to 5' aligned density'''
    
    lengthrange = range(minlength, maxlength + 1)
    
    for length in lengthrange:
        density_plus[length]  = rotate(density_plus[length], -length)
        density_minus[length] = rotate(density_minus[length], length)
        
    return density_plus, density_minus
        
    
#########################################    
##### DATA STROAGE AND MANIPULATION #####
#########################################
    
    
def writebin(datadict, filestring):
    
    f2=open(filestring+"keys","w")
    for chrom in datadict.keys():
        f=open(filestring+chrom,"wb")
        for position in datadict[chrom]:
            f.write(struct.pack("f",position))
        f.close()   
        f2.write(chrom+"\n")
    f2.close()

    
def countstowig(countsfile,filestring):
    import random
    f=open(filestring+".wig","w")
    filestring=filestring.partition("_")[0][-3:]

    random.seed(filestring)
    c1=random.randint(0,255)
    random.seed(filestring+"xxx")
    c2=random.randint(0,255)
    random.seed(filestring+"000")
    c3=random.randint(0,255)

    f.write("track type=wiggle_0 name=tracklabel viewLimits=-5:5 color="+str(c1)+','+str(c2)+','+str(c3)+"\n")
    for chrom in countsfile.keys():
        if chrom[0:3]=='chr':
            f.write("fixedStep  chrom="+chrom+"  start=1  step=1\n")
        else:
            f.write("fixedStep  chrom=\""+chrom+"\"  start=1  step=1\n")

        for position in countsfile[chrom]:
            f.write(str(position)+"\n")
    f.close()
    
def makePickle(data, path_pickle, protocol=pickle.HIGHEST_PROTOCOL):
    f = open(path_pickle, 'w')
    pickle.dump(data, f, protocol=protocol)
    f.close()
    
    
def unPickle(path_pickle):
    #returns the pickled object stored in a pickle file
    f = open(path_pickle, 'r')
    data = pickle.load(f)
    f.close()
    return data


def loadlargePickles(inputs, settings, paths_in, paths_out):
    '''load GFF and density dictionaries into memory to speed up analysis
    WARNING :::: do not run to many samples at once!!!!!'''
    
    files     = inputs['files']
    alignment = settings['alignment']
    minlength = settings['minlength']
    maxlength = settings['maxlength']
    
    lengthindex = range(minlength, maxlength+1)
    
    #load gff
    path_gff_dict = paths_in['path_gff_dict']
    gff_dict = unPickle(path_gff_dict) 
    
    #load density as a larger dict with fname as keys
    path_density = paths_out['path_density']
    plus_dict  = {}
    minus_dict = {}
    
    for fname in files: 
        path_den_plus  = path_density + fname + "/plus_sizesep"
        path_den_minus = path_density + fname + "/minus_sizesep"
        density_plus   = unPickle(path_den_plus)
        density_minus  = unPickle(path_den_minus)
        
        # shift density by readlength for 5' alignment
        if alignment == '5':
            #SD1 == shine dalgarno affinity score
            if fname == 'SD1':    
                for length in lengthindex:
                    density_plus[length] = rotate(density_plus[length], -8) 
                    density_minus[length] = rotate(density_minus[length], 8) 
            else:
                    density_plus, density_minus = alignment_3to5(density_plus, density_minus, minlength, maxlength)

        plus_dict[fname]  = density_plus
        minus_dict[fname] = density_minus

    return gff_dict, plus_dict, minus_dict

    
def shine_dalgarno_affinity(aSD_seq, SD_seq):
    
    '''Calculate diplex energy (Kcal/mol) using RNAstructure:
    Dave Mathews Lab (http://rna.urmc.rochester.edu/RNAstructure.html)'''

    aSD_seq = str(aSD_seq)
    SD_seq  = str(SD_seq)
    
    RNA_prog = RNAstructure(exe_path="/Users/fuad/anaconda/bin/RNAstructure/exe/")
    energy = RNA_prog.DuplexFold(aSD_seq,SD_seq)
    
    return energy


def GC_of_CDS(CDS_seq):
    
    '''Calculate nucleotide composition of each alias in GFF'''
    CDS_length = len(CDS_seq)
    
    G = 0 
    C = 0 
    A = 0 
    T = 0
    
    for nucleotide in CDS_seq:
        if nucleotide == 'G':
            G += 1
        elif nucleotide == 'C':
            C += 1
        elif nucleotide == 'A':
            A += 1
        elif nucleotide == 'T':
            T += 1
    
    G = float(G) / CDS_length * 100
    C = float(C) / CDS_length * 100
    A = float(A) / CDS_length * 100
    T = float(T) / CDS_length * 100
    
    return G, C, A, T
    


def GFF_to_dict(paths_in, gff_settings):
    
    '''Parse gff into dict:
        - feat_of_interest = what to look for in gff (protein_coding, tRNA, rRNA, etc)
        - name_qual        = qualifier for alias/gene name (Name, gene_id)
        - name_qual_alt    = alternative qualifier, if none, set as 'none' 
        - biotype_qual     = qualifier for type of feature (biotype, etc)
        
        These values must correspont to values in the GFF'''
    
    '''Unload gff_settings'''
    
    path_out         = gff_settings['path_out']
    feat_of_interest = gff_settings['feat_of_interest']  #all, protein_coding, tRNA, rRNA
    name_qual        = gff_settings['name_qual']
    name_qual_alt    = gff_settings['name_qual_alt']
    biotype_qual     = gff_settings['biotype_qual']
    aSD_seq          = gff_settings['aSD_seq']
    path_proteins    = paths_in['path_protein']
    
    '''Output path can be defined, or use 0 to set as the annotation file for my main pipeline'''
    
    if path_out == 0:
        path_gff_dict = paths_in['path_gff_dict']
    else:
        path_gff_dict = path_out
    
    '''Parse GFF using BCBio'''
    
    GFFgen = GFF.parse(paths_in['path_gff'])
    chr = GFFgen.next()
    feat_num = 0
    
    '''Define data arrays: will be used as columns for pandas DateFrame'''

    gff_dict   = {}
    aliaslist  = []
    startlist  = []
    stoplist   = []
    seqlist    = []
    strandlist = []
    startcodon = []
    stopcodon  = []
    SDaffinity = []
    G_content  = []
    C_content  = []
    A_content  = []
    T_content  = []

    '''Sift through GFF for relevant information'''
    
    for feature in chr.features:
        
        '''Get feature name'''
        
        if name_qual in feature.qualifiers:
            feat_name = feature.qualifiers[name_qual][0]
        elif name_qual_alt in feature.qualifiers:
            feat_name = feature.qualifiers[name_qual_alt][0]
        else:
            feat_name = 'None'
        
        '''Get feature biotype'''    
        
        if biotype_qual == 'all':
            feat_biotype = 'all' 
        elif biotype_qual == 'protein_coding':   
            
            # for ecoli GFF, no protein_coding qual exists
            # I got protein coding genenames from another database as a .csv:
            
            protein_coding = pd.read_csv(path_proteins)
            proteins = protein_coding.to_dict(orient='list')
            protein_names = proteins['GeneName']
            feat_biotype = 'protein_coding'
        else:
            feat_biotype = 'None'
        
        '''Get start, end, and strand info from GFF'''
        
        start  = chr.features[feat_num].location.start.position 
        end    = chr.features[feat_num].location.end.position
        strand = chr.features[feat_num].strand
                   
        '''Ignore unwanted features'''
       
        if feat_biotype == 'None' or feat_name == 'None':    
            feat_num+=1
            continue
            
            
        if feat_biotype == 'protein_coding':
            if not feat_name in protein_names:    # remove genes not in protein coding list
                feat_num+=1 
                continue
                
        '''Analyze features of interest for feat information'''   
        
        if feat_biotype == feat_of_interest or feat_biotype == 'all':
            
            alias = feat_name
            
            '''Each strand is treated differently, + strand == 1'''
            
            if strand == 1:
                
                '''I save gene sequence + 50 bp from each end:
                makes it easier to analyze start and stop sequence 
                context without using whole genome sequence'''
                
                if start < 50:      # if gene is near the beginning of genome sequence:
                    sequence = 'N' * (50 - start)                  # TB GFF starts at 0, add N * 50
                    sequence = sequence + chr[0:end+50].seq        # gene sequence + 50nt at each end
                else: 
                    sequence = chr[start-50:end+50].seq            # gene sequence + 50nt at each end

                strand_val = '+'
                startcodon_pos = start
                stopcodon_pos  = end-1
                
                if start > 200: 
                    upstream_seq = chr[start-200:start+100].seq
            
            else:
                
                '''For minus strand, 'end' is start codon, 'start' is stop codon
                and sequence is reverse compliment of gene sequence.'''
                
                sequence_rc = chr[start-50:end+50].seq
                sequence = sequence_rc.reverse_complement()
                
                strand_val = '-'
                startcodon_pos = end-1
                stopcodon_pos  = start
                
                if end + 200 > len(chr.seq):
                    upstream_seq = 'none'
                else:
                    upstream_seq_rc = chr[end-100:end+200].seq
                    upstream_seq = upstream_seq_rc.reverse_complement()            
            
            sequence = str(sequence)
            start_codon = sequence[50:53:1]
            stop_codon  = sequence[-53:-50]
            
            '''get sequence from start to stop for GC analysis'''
            
            CDS_seq = sequence[50:-50:1]
            
            G, C, A, T = GC_of_CDS(CDS_seq)
            
            '''Calculate SD affinity'''
            
            SD_seq = sequence[30:50:1]    # analyze 20 nt upstream of start codons
            SD_affinity = shine_dalgarno_affinity(aSD_seq, SD_seq)
            
            '''Append data to lists'''
            
            if alias == 'rpsC':
                print sequence
            
            aliaslist.append(alias)
            seqlist.append(sequence)
            strandlist.append(strand_val)
            startlist.append(startcodon_pos)
            stoplist.append(stopcodon_pos)
            startcodon.append(start_codon)
            stopcodon.append(stop_codon)
            SDaffinity.append(SD_affinity)
            G_content.append(G)
            C_content.append(C)
            A_content.append(A)
            T_content.append(T)
            
            feat_num+=1
            
        else:
            feat_num+=1
            continue
        
    '''Append lists to gff_dict'''

    gff_dict['Alias']    = aliaslist
    gff_dict['Strand']   = strandlist
    gff_dict['Start']    = startlist
    gff_dict['Stop']     = stoplist
    gff_dict['Sequence'] = seqlist
    gff_dict['Start_Codon'] = startcodon
    gff_dict['Stop_Codon']  = stopcodon
    gff_dict['SD_affinity'] = SDaffinity
    gff_dict['G_content'] = G_content
    gff_dict['C_content'] = C_content
    gff_dict['A_content'] = A_content
    gff_dict['T_content'] = T_content

           
    '''Pickle dict for use later'''
    ribo_util.makePickle(gff_dict,path_gff_dict)
    
    '''print dataframe, and save as .csv for use later'''
    ## Print GFF to check 
    gff_df = pd.DataFrame(gff_dict)
    display(gff_df)
    gff_df.to_csv(path_gff_dict + '.csv')
    
    return 
    
    
def heatmapdict_to_df(dictionary, col_name, row_name, val_name):
    '''This is to make dictionary = {index: list_index([values]) into pandas dataframe'''
    
    df = pd.DataFrame(columns = (col_name, row_name, val_name))
    col_index = dictionary.keys()
    
    row_index = [x for x in range(0, len(dictionary[col_index[1]]))]
    rowlist   = []
    collist   = []
    vallist   = []
    
    for column in col_index:
        for row in row_index:
            collist.append(column)
            rowlist.append(row)
            if dictionary[column][row] > 0:
                vallist.append(dictionary[column][row])
            else:
                vallist.append(None)
    
    data = {col_name : collist, 
            row_name : rowlist,
            val_name : vallist
           }
    
    df = pd.DataFrame(data)
    df = df.pivot(col_name, row_name, val_name)
    return df


def dict_to_df(dictionary, col_name, val_name):
    '''This is to make dictionary = {index: list_index(values) into pandas dataframe'''
    
    df = pd.DataFrame(columns = (col_name, val_name))
    col_index = dictionary.keys()
    collist   = []
    vallist   = []
    
    for column in col_index:
        collist.append(column)
        vallist.append(dictionary[column])
        
        
    data = {val_name : vallist}
    df   = pd.DataFrame(data, index = collist)
   
    return df
    

def rotate(l, n):
    return l[-n:] + l[:-n]



    
    
    
    
    
    
    
    
    
    
    
    
    
    