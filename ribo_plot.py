'''
    Script to plot ribosome profiling analysis
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

import pandas as pd
import numpy as np
import math 
from sklearn.linear_model import LinearRegression
import itertools

import plotly 
import plotly.plotly as py
import plotly.graph_objs as go
from plotly.offline import plot
import cPickle as pickle
import seaborn as sns
import matplotlib.pyplot as plt
from IPython.display import display
from scipy import stats

import ribo_util



def size_dist(inputs, paths_in, paths_out):
    
    files = inputs['files']
    path_figure = paths_out['path_figures']

    plot_num = 0
    sns.set_style("white")
    plt.figure(figsize=(5,3))
    
    naming = ''

    for fname in files: 
        
        naming += fname + '_'

        path_analysis = paths_out['path_analysis'] + fname + '/readQC/'
        data = ribo_util.unPickle(path_analysis + 'read_distribution')
        df = ribo_util.dict_to_df(data, 'Length', 'fraction of total')
        plt.plot(df, label = fname)
        plt.title('Size Distribution')
        plt.xlabel("Read Length")
        plt.ylabel("Percent of Reads")
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    plt.savefig(path_figure + 'Comparison/size_dist' + '/sizedist_'+naming+'.pdf', dpi=400, bbox_inches="tight")
    plt.show()

    for fname in files: 
        
        sns.set_style("white")
        plt.figure(figsize=(5,3))


        path_analysis = paths_out['path_analysis'] + fname + '/readQC/'
        data = ribo_util.unPickle(path_analysis + 'read_distribution')
        df = ribo_util.dict_to_df(data, 'Length', 'fraction of total')
        plt.plot(df, label = fname)
        plt.title('Size Distribution')
        plt.xlabel("Read Length")
        plt.ylabel("Percent of Reads")
        plt.legend(loc='upper right')
        
        plt.savefig(path_figure + fname + '/sizedist.pdf', dpi=400, bbox_inches="tight")    
        size_plot_csv = pd.DataFrame(df)
        size_plot_csv.to_csv(path_analysis + 'size_plot_values.csv')
        plt.gcf().clear()
        
def read_comp(inputs, paths_in, paths_out):
    
    files = inputs['files']
    path_figure = paths_out['path_figures']
    
    plot_num = 0
    
    for fname in files: 
        plot_num = 0
        plt.figure(figsize=(20,2.5))
        for nucleotide in ['G', 'C', 'A', 'U']:
            plot_num += 1
            
            path_analysis = paths_out['path_analysis'] + fname + '/readQC/'
            data = ribo_util.unPickle(path_analysis + 'comp_' + nucleotide)
            df = ribo_util.heatmapdict_to_df(data, 'Length', 'Position', 'composition')
            plt.subplot(1,4,plot_num)
            plot = sns.heatmap(df, cmap = "RdBu_r", vmin = 0, vmax = 50)

            plt.setp(plot.get_xticklabels(), visible=False)
            plt.setp(plot.get_xticklabels()[0::5], visible=True)
            plt.setp(plot.get_yticklabels(), visible=False)
            plt.setp(plot.get_yticklabels()[0::4], visible=True)
            plt.title(fname + ' ' + nucleotide)
        plt.savefig(path_figure + fname + '/read_composition.pdf', dpi=400, bbox_inches="tight")  
        plt.show()

def read_comp_context(inputs, paths_in, paths_out):
    
    files = inputs['files']
    path_figure = paths_out['path_figures']
    
    plot_num = 0
    
    for fname in files: 
        plot_num = 0
        plt.figure(figsize=(20,2.5))
        for nucleotide in ['G', 'C', 'A', 'U']:
            plot_num += 1
            
            path_analysis = paths_out['path_analysis'] + fname + '/reads_position/'
            data = ribo_util.unPickle(path_analysis + 'comp_genome_' + nucleotide)
            df = ribo_util.heatmapdict_to_df(data, 'Length', 'Position', 'composition')
            plt.subplot(1,4,plot_num)
            plot = sns.heatmap(df, cmap = "RdBu_r", vmin = 0, vmax = 50)
            
            plot.set_xlim([0,10])

            plt.setp(plot.get_xticklabels(), visible=False)
            plt.setp(plot.get_xticklabels()[0::5], visible=True)
            plt.setp(plot.get_yticklabels(), visible=False)
            plt.setp(plot.get_yticklabels()[0::4], visible=True)
            plt.title(fname + ' ' + nucleotide)
        plt.savefig(path_figure + fname + '/read_comp_context.pdf', dpi=400, bbox_inches="tight")  
        plt.show()


    
'''def read_comp_gff(inputs, paths_in, paths_out):
    
    files = inputs['files']
    
    plot_num = 0
    
    for fname in files: 
        plot_num = 0
        plt.figure(figsize=(20,2.5))
        for nucleotide in ['G', 'C', 'A', 'U']:
            plot_num += 1
            
            path_analysis = paths_out['path_analysis'] + fname + '/'
            data = ribo_util.unPickle(path_analysis + 'comp_genome_' + nucleotide)
            df = ribo_util.heatmapdict_to_df(data, 'Length', 'Position', 'composition')
            plt.subplot(1,4,plot_num)
            plot = sns.heatmap(df, cmap = "RdBu_r", vmin = 0, vmax = 50)

            plt.setp(plot.get_xticklabels(), visible=False)
            plt.setp(plot.get_xticklabels()[0::5], visible=True)
            plt.setp(plot.get_yticklabels(), visible=False)
            plt.setp(plot.get_yticklabels()[0::4], visible=True)
            plt.title(fname + ' ' + nucleotide)
    plt.show()'''
    
    
def plot_alignment_allocation(inputs, paths_in, paths_out):
    files         = inputs['files']
    analysis_path = paths_out['path_analysis_log']
    path_figure   = paths_out['path_figures']

    
    data         = {}
    fnames       = []
    samples      = 0
    samples_list = []
    
    for fname in files: 
    
        fnames.append(fname)
        samples_list.append(samples)
        samples += 1
                
        analysis_log = analysis_path + fname + '/' + fname
        analysis_log = ribo_util.unPickle(analysis_log)

        raw_data = analysis_log['ribo_density']['analysis_breakdown']
        
        total_reads = raw_data.pop('Total Reads', None)
        #total_reads = raw_data.pop('Reads Filtered', None)
                
        for key in raw_data.keys():
            if not key in data:
                data[key] = []
            
            value = raw_data[key]
            data[key].append(value)
            
        legend = raw_data.keys()
        legend = sorted(legend, key=str.lower)

    colors = ['#b2182b','#fddbc7','#e0e0e0','#bababa','#878787','#4d4d4d']
    bottom = [0]*samples
    
    sns.set_style("white")
    plt.figure(figsize=(.8 * samples,4))
 
    for key in legend: 
        i = legend.index(key)
        color = colors[i]
        
        plt.bar(samples_list, data[key], bottom = bottom, color = color, edgecolor='white', width=.6)
    
        bottom = [x + y for x, y in zip(bottom, data[key])]
    plt.xticks(samples_list, fnames, fontweight='bold')
    plt.title('Read Mapping')
    plt.xlabel("Sample")
    plt.ylabel("Reads")
    plt.legend(legend, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    
    plt.savefig(path_figure + fname + '/read_allocation.pdf', dpi=400, bbox_inches="tight")
    plt.show()

    
def plot_avggene(inputs, paths_in, paths_out, settings, settings_plot):
    
    files = inputs['files']
    shift = settings_plot['shift']
    hmmax = settings_plot['HM_max']
    ymax  = settings_plot['ymax']
    
    path_figure = paths_out['path_figures']
    
    minlength      = settings['minlength']
    maxlength      = settings['maxlength']
    length_in_ORF  = settings['length_in_ORF']
    length_out_ORF = settings['length_out_ORF']
    density_type   = settings['density_type'] 
    next_gene      = settings['next_gene']
    equal_weight   = settings['equal_weight']
    threshold      = settings['threshold']
    
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
    
    for fname in files: 
        plot_num = 0
        
        sns.set_style("white")
        sns.set_style("ticks")
        plt.figure(figsize=(20,5))
        path_analysis = paths_out['path_analysis'] + fname + '/avggenes/'
        
        data_start   = ribo_util.unPickle(path_analysis + 'avg_start_'+ name_settings + '_all')
        data_stop    = ribo_util.unPickle(path_analysis + 'avg_stop_' + name_settings + '_all')
        data_startHM = ribo_util.unPickle(path_analysis + 'avg_start_'+ name_settings + '_HM')
        data_stopHM  = ribo_util.unPickle(path_analysis + 'avg_stop_' + name_settings + '_HM')
        
        xmax = len(data_start.keys())
        
        data_start = ribo_util.dict_to_df(data_start, 'Position', 'Reads')
        data_stop  = ribo_util.dict_to_df(data_stop, 'Position', 'Reads')
        data_startHM = ribo_util.heatmapdict_to_df(data_startHM, 'Length', 'Position', 'composition')
        data_stopHM  = ribo_util.heatmapdict_to_df(data_stopHM, 'Length', 'Position', 'composition')
         
        X = range(-50, 151)
        data_startHM = data_startHM.reindex(index=data_startHM.index[::-1])
        data_stopHM  = data_stopHM.reindex(index=data_stopHM.index[::-1])
        
        data_start.to_csv(path_analysis + 'start_all.csv')
        data_stop.to_csv(path_analysis + 'stop_all.csv')
        
        max_start  = data_start["Reads"].max()
        max_stop   = data_stop["Reads"].max()
        
        if ymax == 0:
            if max_start > max_stop:
                ymax = max_start
            else:
                ymax = max_stop
            
        for graph in ['Start', 'Stop']:
            plot_num += 1
            
            if graph == 'Start':
                data = data_start
            elif graph == 'Stop':
                data = data_stop
                
            plt.subplot(2,2,plot_num)
            plt.plot(data, sns.xkcd_rgb["dark grey"], )
            plt.title(fname + ' ' + graph)
            plt.ylabel("Reads")
            plt.ylim(0, ymax)
            plt.xlim(0, xmax)
            sns.despine()
            
        for graph in ['startHM', 'stopHM']:
            plot_num += 1
            if graph == 'startHM':
                dataHM = data_startHM
            elif graph == 'stopHM':
                dataHM = data_stopHM
            
            plt.subplot(2,2,plot_num)
            plot = sns.heatmap(dataHM, cmap = "ocean_r", vmin = 0, vmax = hmmax, cbar=False)

            plt.setp(plot.get_xticklabels(), visible=False)
            plt.setp(plot.get_xticklabels()[0::10], visible=True)
            plt.setp(plot.get_yticklabels(), visible=False)
            plt.setp(plot.get_yticklabels()[0::4], visible=True)
    
        plt.savefig(path_figure + fname + '/avggene_' + name_settings + '.pdf', dpi=400)
    plt.show()
    
    
def plot_frame(inputs, paths_in, paths_out, settings, settings_plot):    
    
    files = inputs['files']
    path_figure = paths_out['path_figures']
    
    minlength = settings['minlength']
    maxlength = settings['maxlength']
    threshold = settings['threshold']
    
    minlength_1      = str(minlength)     +'_'
    maxlength_1      = str(maxlength)     +'_'
    threshold_1      = str(threshold)     +'_'
    
    name_settings = minlength_1+maxlength_1+threshold_1
    
    
    for fname in files:
        sns.set_style("white")
        plt.figure(figsize=(50,3))
        path_analysis = paths_out['path_analysis'] + fname + '/frame/'
        
         
        frame_alias  = ribo_util.unPickle(path_analysis + name_settings + 'frame_alias')        
        alias_F = []
        alias_V = []
        alias_A = []
        alias_R = []
                
        for alias in frame_alias.keys():
            for frame in range(0, 3):
                alias_F.append(frame)
                alias_V.append(frame_alias[alias][frame])
                alias_A.append(alias)
        
        alias_df = pd.DataFrame({'Frame': alias_F, 'Fraction': alias_V, 'Alias': alias_A})
        
        frame_length = ribo_util.unPickle(path_analysis + name_settings + 'frame_length')
        length_F = []
        length_V = []
        length_L = []
        
        for length in frame_length.keys():
            for frame in range(0, 3):
                length_F.append(frame)
                length_V.append(frame_length[length][frame])
                length_L.append(length) 
                
        length_df = pd.DataFrame({'Frame': length_F, 'Fraction': length_V, 'Length': length_L})
        
        frame_genome = ribo_util.unPickle(path_analysis + name_settings + 'frame_genome')
        x_genome = [0, 1, 2]
        
        

        plt.subplot(1,5,1)
        plot1 = sns.barplot(x = x_genome, y = frame_genome, color = ".2")
        plt.ylabel("Fraction")
        plt.xlabel("Frame")
        total = float(len(frame_genome))
        for p in plot1.patches:
            
            height = p.get_height()
            plot1.annotate("%.1f" % p.get_height(), (p.get_x() + p.get_width() / 2., p.get_height()-10),
                 ha='center', va='center', fontsize=10, color='.9', xytext=(0, 20),
                 textcoords='offset points')  

        
        plt.subplot(1,5,3)        
        plot2 = sns.stripplot(x = 'Length', y = 'Fraction', hue = 'Frame', data = length_df,
              jitter=False, size=6, color="0", linewidth=0)
        plot2.legend_.remove() 
        #plot2.axhline(y=45, xmin=0, xmax=1, dashes=[2, 1, 2, 1])

        frame_0 = length_df.loc[length_df.Frame == 0]
     
        x = frame_0.Length.values
        y = frame_0.Fraction.values
        
        plot2.plot(y, sns.xkcd_rgb["pale red"], marker = 'o' )
        
        frame_1 = length_df.loc[length_df.Frame == 1]
        
        x = frame_1.Length.values
        y = frame_1.Fraction.values
        
        plot2.plot(y, c = '0.3', marker = 'o' )

        plt.ylim(0, 100)
        plt.title(fname + ' Frame for each Readlength')
        
        plt.subplot(1,5,5)
        plot3 = sns.boxplot(x = 'Frame', y = 'Fraction', data = alias_df,
                 color="c")
        plt.title(fname + ' Frame for each Alias')
        plot3_1 = sns.stripplot(x = 'Frame', y = 'Fraction', data = alias_df,
              jitter=True, size=3, color=".2", linewidth=0)
    
    plt.savefig(path_figure + fname + '/frame_' + name_settings + '.pdf', dpi=400)
    plt.show()

def plot_frame_comp(inputs, paths_in, paths_out, settings, settings_plot):    
    
    files = inputs['files']
    
    minlength = settings['minlength']
    maxlength = settings['maxlength']
    threshold = settings['threshold']
    
    minlength_1      = str(minlength)     +'_'
    maxlength_1      = str(maxlength)     +'_'
    threshold_1      = str(threshold)     +'_'
    
    name_settings = minlength_1+maxlength_1+threshold_1
    framescore = []
    
    
    for fname in files:
        
        path_analysis = paths_out['path_analysis'] + fname + '/frame/'
        
        frame_genome = ribo_util.unPickle(path_analysis + name_settings + 'frame_genome')
        x_genome = [0, 1, 2]
        frame_score = abs(frame_genome[0] - 40.3)
        frame_score = frame_score / 59.7
        
        framescore.append(frame_score)
        print frame_score 
    print files, framescore
        
    plt.figure(figsize=(30,5))

    plot1 = sns.barplot(x = files, y = framescore, color = ".2")
    plt.ylabel("Libraries")
    plt.xlabel("Frame Score")
    for p in plot1.patches:
            
        height = p.get_height()
        plot1.annotate("%.1f" % p.get_height(), (p.get_x() + p.get_width() / 2., p.get_height()-50),
            ha='center', va='center', fontsize=10, color='0', xytext=(0, 20),
            textcoords='offset points')  
        
        
        
    plt.show()
    
def plot_pausescore(inputs, paths_in, paths_out, settings, settings_plot):
    files      = inputs['files']
    
    aa_codon   = settings_plot['aa_or_codon']
    aminoacids = settings_plot['amino_acid']
    codons     = settings_plot['codon']
    
    ymax_dot   = settings_plot['ymax_dot'] 
    ymax_line  = settings_plot['ymax_line'] 
    vmax_HM    = settings_plot['vmax_HM']   
    
    path_figure = paths_out['path_figures']
    
    aa_plots = len(aminoacids)
    codon_plots = len(codons)
    
    aa_code, codon_code = ribo_util.get_genetic_code()
    
    for fname in files:
        
        minlength       = settings['minlength']
        maxlength       = settings['maxlength']
        plot_upstream   = settings['plot_upstream'] / 3 * 3        #change window to interval of 3
        plot_downstream = settings['plot_downstream'] / 3 * 3
        start_trim      = settings['start_trim'] / 3 * 3
        stop_trim       = settings['stop_trim'] / 3 * 3  
        frameshift      = settings['frameshift']     
        A_site          = settings['A_site shift']   

        
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
    
        path_pausescore = paths_out['path_analysis'] + fname + '/pause_score/'
        
        
        ### For aa_analysis ###
        if aa_codon == 'aa':
            
            aa_score = ribo_util.unPickle(path_pausescore + 'aa_scores' + name_settings) 
            aa_HM    = ribo_util.unPickle(path_pausescore + 'aa_HM_data'+ name_settings ) 
            aa_plot  = ribo_util.unPickle(path_pausescore + 'aa_plot_data'+ name_settings ) 

            aa_df = pd.DataFrame(aa_score)
            aa_df = aa_df.sort_values(by=['Amino Acid'])

            sns.set_style("white")
            sns.set_context("talk")
            plt.figure(figsize=(8 + 4 * aa_plots,5))
            plt.subplot2grid((2,2 + aa_plots), (0,0), rowspan=2, colspan=2)
            plot = sns.stripplot(x="Amino Acid", y="A_site", data=aa_df, size = 13, color= 'black')
            plot = sns.stripplot(x="Amino Acid", y="P_site", data=aa_df, size = 6, color = 'black')
            plot = sns.stripplot(x="Amino Acid", y="E_site", data=aa_df, size = 6, color = 'grey')
            #sns.despine(offset=5, trim = True )
            plt.ylim(0, ymax_dot)
            plot.axhline(y=1, xmin=0, xmax=1, dashes=[2, 2, 2, 2], color = 'grey')

            plt.title(fname + ' Amino Acid Pause Scores')
            plt.xlabel("Amino Acid")
            plt.ylabel("Pause Score")


            aa_HM_dict    = {}
            aa_plot_dict  = {}

            for aa in aminoacids:
                df_HM   = pd.DataFrame(aa_HM[aa])
                df_plot = pd.DataFrame(aa_plot[aa])

                aa_HM_dict[aa]   = df_HM
                aa_plot_dict[aa] = df_plot

            xlim_lower = aa_plot_dict[aa].index[0]
            xlim_upper = aa_plot_dict[aa].index[-1]

            sns.set_style("white")
            sns.set_style("ticks")

            plotnum = 0  
            for aa in aminoacids:

                plt.subplot2grid((2,2 + aa_plots), (0,2 + plotnum), rowspan=1, colspan=1)

                plot_2 = plt.plot(aa_plot_dict[aa], sns.xkcd_rgb["dark grey"])
                plt.title(aa + ' Average Plot')
                plt.ylabel("Pause Score")
                plt.ylim(0, ymax_line)
                plt.xlim(xlim_lower, xlim_upper)
                sns.despine(offset=5)


                plt.subplot2grid((2,2 + aa_plots), (1,2 + plotnum), rowspan=1, colspan=1)
                plot_3 = sns.heatmap(aa_HM_dict[aa], cmap = "ocean_r", vmin = 0, vmax = vmax_HM,
                                     cbar=False, xticklabels=15, yticklabels=6)
                sns.despine(offset=5)

                plotnum +=1

            plt.tight_layout()
            plt.savefig(path_figure + fname + '/aa_pausescore' + name_settings + '_aa_pause_scores.pdf', dpi=400)
            plt.show()
            
            aa_plot_csv = pd.DataFrame(aa_plot)
            aa_plot_csv.to_csv(path_pausescore + 'aa_plot_values.csv')

        ### For codon_analysis ###
        if aa_codon == 'codon':
            
            codon_score = ribo_util.unPickle(path_pausescore + 'codon_scores' + name_settings ) 
            codon_HM    = ribo_util.unPickle(path_pausescore + 'codon_HM_data' + name_settings ) 
            codon_plot  = ribo_util.unPickle(path_pausescore + 'codon_plot_data' + name_settings ) 
            
            aa_list = []
            for codon in codon_score['Codon']:
                aa = codon_code[codon]
                aa_list.append(aa)
            
            codon_score['Amino_Acid'] = aa_list    
                
            codon_df = pd.DataFrame(codon_score)
            codon_df = codon_df.sort_values(by=['Amino_Acid', 'Codon'])

            sns.set_style("white")
            sns.set_context("talk")
            plt.figure(figsize=(28 + 4 * codon_plots,5))
            plt.subplot2grid((2,7 + codon_plots), (0,0), rowspan=2, colspan=7)
            
            plot = sns.stripplot(x="Codon", y="A_site", data=codon_df, size = 10, color = 'black')
            plot = sns.stripplot(x="Codon", y="P_site", data=codon_df, size = 2, color = 'black')
            plot = sns.stripplot(x="Codon", y="E_site", data=codon_df, size = 2, color = 'grey')
            #sns.despine(offset=5, trim = True )
            plt.ylim(0, ymax_dot)
            plot.axhline(y=1, xmin=0, xmax=1, dashes=[2, 2, 2, 2], color = 'grey') 
            
            plot.axvline(x=3.5, ymin=0, ymax=ymax_dot, color = 'grey', lw = .8)
            plt.text(1, 4,'A', fontsize=30) #add text
            plot.axvline(x=5.5, ymin=0, ymax=ymax_dot, color = 'grey', lw = .8)
            plt.text(4, 4,'C', fontsize=30) #add text
            plot.axvline(x=7.5, ymin=0, ymax=ymax_dot, color = 'grey', lw = .8)
            plt.text(6, 4,'D', fontsize=30) #add text
            plot.axvline(x=9.5, ymin=0, ymax=ymax_dot, color = 'grey', lw = .8)
            plt.text(8, 4,'E', fontsize=30) #add text
            plot.axvline(x=11.5, ymin=0, ymax=ymax_dot, color = 'grey', lw = .8)
            plt.text(10, 4,'F', fontsize=30) #add text
            plot.axvline(x=15.5, ymin=0, ymax=ymax_dot, color = 'grey', lw = .8)
            plt.text(13, 4,'G', fontsize=30) #add text
            plot.axvline(x=17.5, ymin=0, ymax=ymax_dot, color = 'grey', lw = .8)
            plt.text(16, 4,'H', fontsize=30) #add text
            plot.axvline(x=20.5, ymin=0, ymax=ymax_dot, color = 'grey', lw = .8)
            plt.text(18.8, 4,'I', fontsize=30) #add text
            plot.axvline(x=22.5, ymin=0, ymax=ymax_dot, color = 'grey', lw = .8)
            plt.text(21, 4,'K', fontsize=30) #add text
            plot.axvline(x=28.5, ymin=0, ymax=ymax_dot, color = 'grey', lw = .8)
            plt.text(25, 4,'L', fontsize=30) #add text
            plot.axvline(x=29.5, ymin=0, ymax=ymax_dot, color = 'grey', lw = .8)
            plt.text(28.6, 4,'M', fontsize=30) #add text
            plot.axvline(x=31.5, ymin=0, ymax=ymax_dot, color = 'grey', lw = .8)
            plt.text(30, 4,'N', fontsize=30) #add text
            plot.axvline(x=35.5, ymin=0, ymax=ymax_dot, color = 'grey', lw = .8)
            plt.text(33, 4,'P', fontsize=30) #add text
            plot.axvline(x=37.5, ymin=0, ymax=ymax_dot, color = 'grey', lw = .8)
            plt.text(36, 4,'Q', fontsize=30) #add text
            plot.axvline(x=43.5, ymin=0, ymax=ymax_dot, color = 'grey', lw = .8)
            plt.text(40, 4,'R', fontsize=30) #add text
            plot.axvline(x=49.5, ymin=0, ymax=ymax_dot, color = 'grey', lw = .8)
            plt.text(46, 4,'S', fontsize=30) #add text
            plot.axvline(x=53.5, ymin=0, ymax=ymax_dot, color = 'grey', lw = .8)
            plt.text(51, 4,'T', fontsize=30) #add text
            plot.axvline(x=57.5, ymin=0, ymax=ymax_dot, color = 'grey', lw = .8)
            plt.text(55, 4,'V', fontsize=30) #add text
            plot.axvline(x=60.5, ymin=0, ymax=ymax_dot, color = 'grey', lw = .8)
            plt.text(58.4, 4,'W', fontsize=30) #add text
            plt.text(61.5, 4,'_', fontsize=30) #add text'''
            
            
            
            
            plt.title(fname + ' Codon Pause Scores')
            plt.xlabel("Codon")
            plt.ylabel("Pause Score")


            codon_HM_dict    = {}
            codon_plot_dict  = {}

            for codon in codons:
                df_HM   = pd.DataFrame(codon_HM[codon])
                df_plot = pd.DataFrame(codon_plot[codon])

                codon_HM_dict[codon]   = df_HM
                codon_plot_dict[codon] = df_plot

            xlim_lower = codon_plot_dict[codon].index[0]
            xlim_upper = codon_plot_dict[codon].index[-1]

            sns.set_style("white")
            sns.set_style("ticks")

            plotnum = 0  
            for codon in codons:

                plt.subplot2grid((2,7 + codon_plots), (0,7 + plotnum), rowspan=1, colspan=1)

                plot_2 = plt.plot(codon_plot_dict[codon], sns.xkcd_rgb["dark grey"])
                plt.title(codon + ' Average Plot')
                plt.ylabel("Pause Score")
                plt.ylim(0, ymax_line)
                plt.xlim(xlim_lower, xlim_upper)
                sns.despine(offset=5)


                plt.subplot2grid((2,7 + codon_plots), (1,7 + plotnum), rowspan=1, colspan=1)
                plot_3 = sns.heatmap(codon_HM_dict[codon], cmap = "ocean_r", vmin = 0, vmax = vmax_HM,
                                     cbar=False, xticklabels=15, yticklabels=6)
                sns.despine(offset=5)

                plotnum +=1



            plt.tight_layout()
            plt.savefig(path_figure + fname + '/codon_pausescore' + name_settings + '.pdf', dpi=400)
            plt.show()

def plot_asymmetry(inputs, paths_in, paths_out, settings):    
    
    files = inputs['files']
    path_figure = paths_out['path_figures']
    
    minlength    = settings['minlength']
    maxlength    = settings['maxlength']
    
    minlength_1       = str(minlength)      +'_'
    maxlength_1       = str(maxlength)      +'_'
    
    name_settings = minlength_1+maxlength_1
    
    for fname in files:

        path_analysis = paths_out['path_analysis'] + fname + '/'
        
        asymmetry_dict  = ribo_util.unPickle(path_analysis + 'asymmetry') 
        df = pd.DataFrame(asymmetry_dict)
        #display(df)
        sub   = df.loc[df.Subgroup == "Subgroup"]
        other = df.loc[df.Subgroup == "Other"]
        
         
        
        plot2 = sns.jointplot(x="Score", y="GeneLength", data=df, size = 4, 
                    color="black", marker='.', alpha=0, xlim = (-3, 3))
        
        plot2.x = other.Score
        plot2.y = other.GeneLength
        
        plot2.plot_joint(plt.scatter, marker='.', c='.6')
        plot2.plot_joint(sns.kdeplot, n_levels=6, cmap="Greys")
        
        plot2.x = sub.Score
        plot2.y = sub.GeneLength

        plot2.plot_joint(plt.scatter, marker='.', c='orangered', alpha=.6)
        plot2.plot_joint(sns.kdeplot, n_levels=4, cmap="Reds")
        
        plt.savefig(path_figure + fname + '/asymmetry' + name_settings + '.pdf', dpi=400, bbox_inches="tight") 
        plt.show()
    
def plot_asymmetry_comp(inputs, paths_in, paths_out, settings):
    #sns.boxplot(data=iris, orient="h", palette="Set2")
    files = inputs['files']
    path_figure = paths_out['path_figures'] + 'Comparison/'
    
    minlength    = settings['minlength']
    maxlength    = settings['maxlength']
    
    minlength_1       = str(minlength)      +'_'
    maxlength_1       = str(maxlength)      +'_'
    
    name_settings = minlength_1+maxlength_1
    
    data  = {}
    score_list = []
    fname_list = []
    
    plot_count = 0
    fnames = ''
    cm_treated = ['CW63', 'OL38', 'OL34', 'OL35', 'OL18', 'OL19', 'OL20', 'OL21', 'OL1', 'OL2']
    cm_untreated = ['OL' + str(i) for i in range(0, 40)]
    
    treated_mean = []
    untreated_mean = []
    
    for fname in files:
        
            
        fnames += fname
        plot_count += 1

        path_analysis = paths_out['path_analysis'] + fname + '/'
        asymmetry_dict  = ribo_util.unPickle(path_analysis + 'asymmetry') 
        df = pd.DataFrame(asymmetry_dict)
        
        compare_path = paths_out['path_analysis'] + 'OL17' + '/'
        compare_dict  = ribo_util.unPickle(compare_path + 'asymmetry') 
        compare_df = pd.DataFrame(compare_dict)
        
        score = df['Score'].values.tolist()
        compare_score = compare_df['Score'].values.tolist()
        
        ttest = stats.ttest_ind(score, compare_score )
        ttest2 = stats.ttest_1samp(score, 0 )
        
        mean = np.mean(score) 
        
        if fname in cm_treated:
            treated_mean.append(mean)
        elif fname in cm_untreated:
            untreated_mean.append(mean)
                
        for value in score: 
            fname_list.append(fname)
            score_list.append(value)
    
    data['file'] = fname_list
    data['score'] = score_list
    data = pd.DataFrame(data)
    
    
    sns.set_style("white")
    plot = plt.figure(figsize=(3,12))
    plot = sns.boxplot(x="score", y='file', data=data, orient="h", color = "grey", fliersize = 0)
    plot.axvline(linewidth=2, color='r')
    plot.set(xlim=(-3, 3))
    sns.despine(offset=5)
    plt.savefig(path_figure + 'asymmetry/' + fnames + name_settings + '.pdf', dpi=400, bbox_inches="tight") 

    plt.show()
    
    plot2 = sns.distplot(treated_mean, bins=5)
    plot2 = sns.distplot(untreated_mean, bins=5)
    plt.show()
    
    ttest = stats.ttest_ind(treated_mean, untreated_mean )
    
    print ttest
    print 'number of genes in treated = ' + str(len(treated_mean))
    print 'number of genes in untreated = ' + str(len(untreated_mean))
    print
    print 'treated = ' , treated_mean
    print
    print 'untreated = ' , untreated_mean


def plot_genelist(fnamelist, paths_in, paths_out, settings):
    '''Takes 2 files and compares gene RPKM values'''
    path_figure = paths_out['path_figures'] + 'Comparison/'

    '''Load Settings'''
    file1 = fnamelist[0]
    file2 = fnamelist[1]
    foldchange      = settings['foldchange']
    highlight_genes = settings['interesting']
    data_transform  = settings['data_transform'] # has to be 'None', 'log10', or 'log2' 
    highlight_genes = settings['annotation']
    
    '''Load Data'''
    
    #genelist file from analysis
    path_analysis1 = paths_out['path_analysis'] + fnamelist[0] + '/'
    path_analysis2 = paths_out['path_analysis'] + fnamelist[1] + '/'
   
    genelist_f1 = ribo_util.unPickle(path_analysis1 + 'genelist') 
    genelist_f2 = ribo_util.unPickle(path_analysis2 + 'genelist') 

    #define alias and RPKM as arrays
    df1           = pd.DataFrame(genelist_f1)
    alias         = df1.Alias
    RPKM_f1       = df1.RPKM
    annotation_f1 = df1.Annotation
    
    df2           = pd.DataFrame(genelist_f2)
    RPKM_f2       = df2.RPKM
    annotation_f2 = df2.Annotation
    
    '''data arrays to plot'''
    xydict = {}
    x_list = []
    y_list = []
    folddict = {}
    x_foldchange  = []
    y_foldchange  = []
    interestingdict = {}
    x_interesting = []
    y_interesting = []
    annotateddict = {}
    x_annotated = []
    y_annotated = []
        
    '''iterate through gene data for manipulation (log transform and fold change calc)'''    
    
    for alias, xval, yval, annotation in itertools.izip(alias.values, RPKM_f1.values, RPKM_f2.values, annotation_f1.values):
        
        # remove genes with RPKM = 0
        if xval == 0 or yval == 0:
            continue
        
        # calculate foldchange prior to data transformation
        foldchange_val = (yval / xval)
       
        # transform data
        if data_transform == 'None':
            xval = xval
            yval = yval
        if data_transform == 'log2':
            xval = math.log(xval, 2)
            yval = math.log(yval, 2)
        if data_transform == 'log10':
            xval = math.log10(xval)
            yval = math.log10(yval)
            
        # append data to lists  
        
        x_list.append(xval)
        y_list.append(yval)
        
        if foldchange_val > foldchange or foldchange_val < (1 / foldchange):
            x_foldchange.append(xval)
            y_foldchange.append(yval)
        
        if alias in settings['interesting']:
            x_interesting.append(xval)
            y_interesting.append(yval)
            
            print alias, xval, yval
                        
        if annotation == highlight_genes:
            x_annotated.append(xval)
            y_annotated.append(yval)
            
            print alias, xval, yval
    
    # convert lists to pd.dataframe
    
    xydict[file1] = x_list
    xydict[file2] = y_list
    xy_df = pd.DataFrame(xydict)
    
    folddict[file1] = x_foldchange
    folddict[file2] = y_foldchange
    fold_df = pd.DataFrame(folddict)
    
    interestingdict[file1] = x_interesting
    interestingdict[file2] = y_interesting
    interesting_df = pd.DataFrame(interestingdict)
    
    annotateddict[file1] = x_annotated
    annotateddict[file2] = y_annotated
    annotated_df = pd.DataFrame(annotateddict)
    
    display(annotated_df)
    
    sns.set_style("white")
    plt.figure(figsize=(6,6))

    plot1 = sns.regplot(x = file1, y = file2, data = xy_df, 
                            color=".5", marker='.', fit_reg = False)
    
    #plot1.set_axis_labels(fnamelist[0] + ' log2(RPKM)', fnamelist[1] + ' log2(RPKM)')
    
    plot1.x = fold_df[file1]
    plot1.y = fold_df[file2]
        
    #plot1.plot_joint(plt.scatter, marker='.', c='r', alpha = .4)
    
    x = interesting_df[file1]
    y = interesting_df[file2]
    
    plot1 = sns.regplot(x = x, y = y, color="blue", marker='.', fit_reg = False)
    
    x = annotated_df[file1]
    y = annotated_df[file2]
    
    plot1 = sns.regplot(x = x, y = y, color="red", marker='.', fit_reg = False)
    
    x = xy_df[file1]
    y = xy_df[file1]
    
    plot1 = sns.regplot(x = x, y = y, color="black", marker="None")
    
    plot2 = sns.jointplot(x = file1, y = file2, data = xy_df, 
                            color=".5", marker='.', kind="reg")
    plot2 = plot2.annotate(stats.pearsonr)

    plt.savefig(path_figure + file1 + '_' + file2 + '.pdf', dpi=400, bbox_inches="tight") 
    

    plt.show()
    
    
def plot_UTRvsCDS(inputs, paths_in, paths_out, settings, settings_plot):
    
    files = inputs['files']
    path_figure = paths_out['path_figures']
    
    minlength      = settings['minlength']
    maxlength      = settings['maxlength']
    length_in_ORF  = settings['length_in_ORF']
    length_out_ORF = settings['length_out_ORF']
    density_type   = settings['density_type'] 
    next_gene      = settings['next_gene']
    equal_weight   = settings['equal_weight']
    threshold      = settings['threshold']
    
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
    
    all_data = {}
    sns.set_style("white")
    plot = plt.figure(figsize=(5,5))
    
    for fname in files: 
        
        path_analysis = paths_out['path_analysis'] + fname + '/avggenes/'
        
        data = ribo_util.unPickle(path_analysis + 'avg_start_'+ name_settings + '_UTR_CDS')
        #all_data[fname] = data.values()
        
        total_values = len(data.values())
        sum_values   = sum(data.values())
        
        print fname, total_values
        
        df_sort = [ i / sum_values for i in data.values()]              
        df = pd.Series(data.values())
        df = df.sort_values(ascending=True)
        
        plot = plt.hist(df, total_values, normed=1, histtype='step',cumulative=True, label = fname)
        plt.xlim(-1, 2)
        
        plt.legend(loc='right')
        
    plt.show()
    
    
    '''all_data = pd.DataFrame(all_data)
    plot = plt.figure(figsize=(5,5))
    plot = sns.kdeplot(all_data, cumulative=True)'''
        
    
    
    
    
def plot_TEvsSD(fnamelist, paths_in, paths_out, settings):
    '''Takes 2 files and compares translation efficiency 
    to shine dalgarno affinity'''
    
    '''Load Settings'''
    
    foldchange      = settings['foldchange']
    highlight_genes = settings['interesting']
    data_transform  = settings['data_transform'] # has to be 'None', 'log10', or 'log2' 
    
    '''Load Data'''
    
    #genelist file from analysis
    path_analysis1 = paths_out['path_analysis'] + fnamelist[0] + '/'
    path_analysis2 = paths_out['path_analysis'] + fnamelist[1] + '/'
   
    genelist_f1 = ribo_util.unPickle(path_analysis1 + 'genelist') 
    genelist_f2 = ribo_util.unPickle(path_analysis2 + 'genelist') 

    #define alias and RPKM as arrays
    df1      = pd.DataFrame(genelist_f1)
    alias    = df1.Alias
    RPKM_f1  = df1.RPKM
    
    df2 = pd.DataFrame(genelist_f2)
    RPKM_f2   = df2.RPKM
    
    '''data arrays to plot'''
    
    x = []
    y = []
    x_foldchange  = []
    y_foldchange  = []
    x_interesting = []
    y_interesting = []
        
    '''iterate through gene data for manipulation (log transform and fold change calc)'''    
    
    for alias, xval, yval in itertools.izip(alias.values, RPKM_f1.values, RPKM_f2.values):
        
        # remove genes with RPKM = 0
        if xval == 0 or yval == 0:
            continue
        
        # calculate foldchange prior to data transformation
        foldchange = abs(yval - xval)
       
        # transform data
        if data_transform == 'None':
            xval = xval
            yval = yval
        if data_transform == 'log2':
            xval = math.log2(xval)
            yval = math.log2(yval)
        if data_transform == 'log10':
            xval = math.log10(xval)
            yval = math.log10(yval)
            
        # append data to lists  
        
        x.append(xval)
        y.append(yval)
        
        if foldchange > settings['foldchange']:
            x_foldchange.append(xval)
            y_foldchange.append(yval)
            
        if alias in settings['interesting']:
            x_interesting.append(xval)
            y_interesting.append(yval)
            
    plot1 = sns.jointplot(x=x, y=y, size = 6, 
                            color=".4", marker='.', alpha = .5)
    
    plot1.set_axis_labels(fnamelist[0] + ' log2(RPKM)', fnamelist[1] + ' log2(RPKM)')
    
    plot1.x = x_foldchange
    plot1.y = y_foldchange
        
    plot1.plot_joint(plt.scatter, marker='.', c='r', alpha = .4)
    
    plot1.x = x_interesting
    plot1.y = y_interesting
    
    plot1.plot_joint(plt.scatter, marker='.', c='b', alpha = .8)

    
    plt.show()
        
    
    
    
    
    
 
    
def plot_avggene_end(inputs, paths_in, paths_out, settings):
    
    files = inputs['files']
    shift = settings['shift']
    hmmax = settings['HM_max']
    
    for fname in files: 
        plot_num = 0
        plt.figure(figsize=(20,5))
        path_analysis = paths_out['path_analysis'] + fname + '/'
        
        data_start = ribo_util.unPickle(path_analysis + 'avg_start_all_end')
        data_stop  = ribo_util.unPickle(path_analysis + 'avg_stop_all_end')
        data_startHM = ribo_util.unPickle(path_analysis + 'avg_start_HM_end')
        data_stopHM  = ribo_util.unPickle(path_analysis + 'avg_stop_HM_end')
        
        xmax = len(data_start.keys())
        
        data_start = ribo_util.dict_to_df(data_start, 'Position', 'Reads')
        data_stop  = ribo_util.dict_to_df(data_stop, 'Position', 'Reads')
        data_startHM = ribo_util.heatmapdict_to_df(data_startHM, 'Length', 'Position', 'composition')
        data_stopHM  = ribo_util.heatmapdict_to_df(data_stopHM, 'Length', 'Position', 'composition')
        
        data_start.to_csv(path_analysis + 'start_all.csv')
        data_stop.to_csv(path_analysis + 'stop_all.csv')
        
        max_start  = data_start["Reads"].max()
        max_stop   = data_stop["Reads"].max()
        
        if max_start > max_stop:
            ymax = max_start
        else:
            ymax = max_stop
            
        for graph in ['Start', 'Stop']:
            plot_num += 1
            
            if graph == 'Start':
                data = data_start
            elif graph == 'Stop':
                data = data_stop
                
            plt.subplot(2,2,plot_num)
            plt.plot(data, sns.xkcd_rgb["dark grey"], )
            plt.title(fname + ' ' + graph)
            plt.ylabel("Reads")
            plt.ylim(0, ymax)
            plt.xlim(0, xmax)
            sns.despine()
        for graph in ['startHM', 'stopHM']:
            plot_num += 1
            if graph == 'startHM':
                dataHM = data_startHM
            elif graph == 'stopHM':
                dataHM = data_stopHM
            
            plt.subplot(2,2,plot_num)
            plot = sns.heatmap(dataHM, cmap = "ocean_r", vmin = 0, vmax = hmmax, cbar=False)

            plt.setp(plot.get_xticklabels(), visible=False)
            plt.setp(plot.get_xticklabels()[0::10], visible=True)
            plt.setp(plot.get_yticklabels(), visible=False)
            plt.setp(plot.get_yticklabels()[0::4], visible=True)
            
    plt.show()
    
    
    
def plot_pausescore_APE_heatmap(inputs, paths_in, paths_out, settings, settings_plot):
    files      = inputs['files']
    library_id = pd.read_csv(paths_in['files'])

    
    center_HM = settings_plot['center_HM'] 
    vmax_HM   = settings_plot['vmax_HM']   
    
    
    namelist = []
    library  = []
    Amino_acid = []
    A_score = []
    P_score = []
    E_score = []
    
    A_dict = {}
    P_dict = {}
    E_dict = {}
    all_dict = {}
    plots = 0 
    
    
    for fname in library_id.Library:
        lib        = library_id.loc[library_id.Library == fname]
        lib_index  = lib.Name.index[0]
        name       = lib.Name.loc[lib_index]
        sort       = lib.Sort.loc[lib_index]
        
        minlength       = settings['minlength']
        maxlength       = settings['maxlength']
        plot_upstream   = settings['plot_upstream'] / 3 * 3        #change window to interval of 3
        plot_downstream = settings['plot_downstream'] / 3 * 3
        start_trim      = settings['start_trim'] / 3 * 3
        stop_trim       = settings['stop_trim'] / 3 * 3  
        frameshift      = settings['frameshift']     

        
        minlength_1       = str(minlength)      +'_'
        maxlength_1       = str(maxlength)      +'_'
        plot_upstream_1   = str(plot_upstream)  +'_'
        plot_downstream_1 = str(plot_downstream)+'_'
        start_trim_1      = str(start_trim)     +'_'
        stop_trim_1       = str(stop_trim)      +'_'
        frameshift_1      = str(frameshift)     +'_'

        name_settings = plot_upstream_1+plot_downstream_1+start_trim_1+stop_trim_1+minlength_1+maxlength_1+frameshift_1
        
        path_pausescore = paths_out['path_analysis'] + fname + '/pause_score/'
        
        aa_score = ribo_util.unPickle(path_pausescore + name_settings + 'aa_scores') 
        aa_HM    = ribo_util.unPickle(path_pausescore + name_settings + 'aa_HM_data') 
        aa_plot  = ribo_util.unPickle(path_pausescore + name_settings + 'aa_plot_data')
        
        amino_acids = aa_score['Amino Acid']
        a_site      = aa_score['A_site']
        p_site      = aa_score['P_site']
        e_site      = aa_score['E_site']
        
        plots += 1
        for aa, ascore, pscore, escore in itertools.izip(amino_acids, a_site ,p_site, e_site): 
            
            if aa == '_':
                continue 
            namelist.append(fname)
            library.append(sort)
            Amino_acid.append(aa)
            A_score.append(ascore)
            P_score.append(pscore)
            E_score.append(escore)

            
    A_dict['Library']    = library
    A_dict['Amino Acid'] = Amino_acid
    A_dict['A site']     = A_score
    
    P_dict['Library']    = library
    P_dict['Amino Acid'] = Amino_acid
    P_dict['P site']     = P_score
    
    E_dict['Library']    = library
    E_dict['Amino Acid'] = Amino_acid
    E_dict['E site']     = E_score
    
    all_dict['Library']    = namelist
    all_dict['Amino Acid'] = Amino_acid
    all_dict['A site']     = A_score
    all_dict['P site']     = P_score
    all_dict['E site']     = E_score
    
    A_df = pd.DataFrame(A_dict)
    P_df = pd.DataFrame(P_dict)
    E_df = pd.DataFrame(E_dict)
    all_df = pd.DataFrame(all_dict)

    A_df = A_df.pivot('Library', 'Amino Acid', 'A site')
    P_df = P_df.pivot('Library', 'Amino Acid', 'P site')
    E_df = E_df.pivot('Library', 'Amino Acid', 'E site')
    all_df = all_df[['Library', 'Amino Acid', 'A site', 'P site', 'E site']]    
       
    sns.set_style("white")
    sns.set_context("talk")
    plt.figure(figsize=(20,.8+.4* plots))
    
    
    plt.subplot(1, 3, 1)
    plot_A = sns.heatmap(A_df, cmap = "RdBu_r", vmax = vmax_HM, center = center_HM,
                             cbar=False,)
    plt.title('A-Site Pause Scores')
    
    plt.subplot(1, 3, 2)
    plot_P = sns.heatmap(P_df, cmap = "RdBu_r", vmax = vmax_HM, center = center_HM,
                             cbar=False, )
    plt.title('P-Site Pause Scores')
    
    plt.subplot(1, 3, 3)
    plot_E = sns.heatmap(E_df, cmap = "RdBu_r", vmax = vmax_HM, center = center_HM,
                             cbar=True, )
    plt.title('E-Site Pause Scores')
    
    plt.tight_layout()  
    
    outpath = paths_out['path_figures']
    plt.savefig(outpath + '/Comparison/' + 'pause_scores.png', dpi=400)
    plt.show()
    
    
    #display(all_df)
    #all_df.to_csv(inpath + '2codon_scores.csv')
    
def plot_pausescore_downstream(inputs, paths_in, paths_out, settings, settings_plot):
    files      = inputs['files']
    library_id = pd.read_csv(paths_in['files'])

    aminoacids = settings_plot['amino_acid']

    center_HM = settings_plot['center_HM'] 
    vmax_HM   = settings_plot['vmax_HM']   
      
    
    datadict = {}
    
    lib_names  = []
    plot_data  = []
    amino_acid = []
    position   = []
    
    plots = 0 
    for fname in library_id.Library:
        lib        = library_id.loc[library_id.Library == fname]
        lib_index  = lib.Name.index[0]
        name       = lib.Name.loc[lib_index]
        sort       = lib.Sort.loc[lib_index]
        
        minlength       = settings['minlength']
        maxlength       = settings['maxlength']
        plot_upstream   = settings['plot_upstream_wave'] / 3 * 3        #change window to interval of 3
        plot_downstream = settings['plot_downstream_wave'] / 3 * 3
        start_trim      = settings['start_trim'] / 3 * 3
        stop_trim       = settings['stop_trim'] / 3 * 3  
        frameshift      = settings['frameshift']   
        next_codon      = settings['next_codon']

        
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
        
        path_pausescore = paths_out['path_analysis'] + fname + '/pause_score/waves/'
        
  
        aa_plot  = ribo_util.unPickle(path_pausescore + name_settings + 'aa_plot_data')
        df_plot = pd.DataFrame(aa_plot)
        
        for aa in aminoacids:
            
            data = df_plot[aa].tolist()
            data_length = len(data)
            data_index  = 0
            
            
            for index in range(0, data_length/3):
                codon_value = data[data_index-1: data_index+2]
                codon_value = sum(codon_value) / 3

                lib_names.append(sort)
                plot_data.append(codon_value)
                amino_acid.append(aa)
                position.append(index)
                data_index += 3
        
        plots += 1
            
        
    datadict['Sort'] = lib_names
    datadict['Data'] = plot_data
    datadict['AA']   = amino_acid
    datadict['Position']   = position    
        
    data_df = pd.DataFrame(datadict)      
        
    for aa in aminoacids: 
        
        aa_data = data_df.loc[data_df.AA == aa]
        aa_data = aa_data.drop('AA', axis=1)
        
        aa_data = aa_data.pivot('Sort', 'Position', 'Data')
    
        sns.set_style("white")
        sns.set_context("talk")
        plt.figure(figsize=(20,.8+.4* plots))
        
        vmax   = vmax_HM
        center = center_HM
        vmin   = center_HM - (vmax_HM-center)
        
        plot_A = sns.heatmap(aa_data, cmap = "RdBu_r", vmax = vmax, center = center, vmin = vmin,
                                 cbar=True,)
        plt.title(aa)

        
        outpath = paths_out['path_figures']
        plt.savefig(outpath + '/Comparison/' + aa + 'pause_wave.png', dpi=400)
        plt.show()
    
    #all_df.to_csv(inpath + 'codon_scores.csv')

        