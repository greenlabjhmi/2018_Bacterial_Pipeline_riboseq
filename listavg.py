# This file contains code for creating average plots and creating gene lists. 
# The workflow function automates a lot of it, but still can be tweaked depending on needs.
from Bio import SeqIO
from Bio.Seq import Seq
import csv
#import searchcode

# This is the workflow that does the quantitation of the genes for ORF, 5'UTR, and 3'UTR.
# bp5 and bp3 are extra distance on the 5' and 3' ends of the "gene" that will be included in quantitation.
# Gene is defined as ORF is the UTR GFFs are entered as -1s. Otherwise, it will include the UTRs, if available. 
# This workflow will go through the genome and generate a gene list. Filtermodules can be called to select for particular genes.
# filebase is the name that will be used in file creation along with the file path.
# counts_filestring is path and rootname for density lists (plus and minus).
# shift is amount ribosomes are shifted for list creation.
# filter is T or F, do you want to filter data, done in function and filtermodule
# GFFgen_filename is a string of the filename with full path
# thresh signifies minimal rpkm needed in coding region for gene to be in list.

def totalquant_wf(filebase,counts_filestring,bp5,bp3,ignoreutr5,ignoreutr3,shift,filter,
				thresh,GFFgen_filename,utrgfffilename,utr5gfffilename,totalreads):

# 	comments="This is a test run of totalquant_wf.\nIt makes a list (see csv file).\n"
# 	comments+="shift is amount reads are shifted for calculation of gene counts in the list.\n"
# 	comments+="bp5 and bp3 are the extra distance on 5 and 3' ends for average. It is the entire UTR if no UTR info available for that gene. -1 for UTR means that all UTRs will be bp distance.\n"
# 	comments+="thresh signifies minimal rpkm needed in coding region for gene to be in list.\n"
# 	# Write output file of comments.
# 	fc=open(filebase+"_output.txt","w")
# 	fc.write(comments)
# 	fc.write("\n")
# 	fc.write("Totalgene_wf was called with parameters:\n")
# 	fc.write("counts_filestring="+str(counts_filestring)+"\n")
# 	fc.write("bp5,bp3="+str(bp5)+", "+str(bp3)+"\n")
# 	fc.write("shift="+str(shift)+"\n")
# 	fc.write("filter="+str(filter)+"\n")
# 	fc.write("thresh="+str(thresh)+"\n")
# 	fc.write("ignoreutr5="+str(ignoreutr5)+"\n")
# 	fc.write("ignoreutr3="+str(ignoreutr3)+"\n")
# 	fc.close()

	GFFgen=GFF.parse(GFFgen_filename)
	GFFlist=seqtools.makeGFFlist(GFFgen)

	counts0=builddense.readcountsf(counts_filestring+"_plus_")		
	counts1=builddense.readcountsf(counts_filestring+"_minus_")
	counts=[counts0,counts1]
	
	if ignoreutr3 != 1:
		utrgffgen = GFF.parse(utrgfffilename)
		utrtable = genometools.makeutrtable(utrgffgen)
	else:
		utrtable = {}
	
	if ignoreutr5 != 1:
		utrgffgen = GFF.parse(utr5gfffilename)
		utrtable2 = genometools.makeutrtable(utrgffgen)
	else:
		utrtable2 = {}
				
	mgl=makegenelist(counts,GFFlist,utrtable,utrtable2,bp5,bp3,
					ignoreutr5,ignoreutr3,shift,filter,thresh,totalreads)
	writedicttoexcel(mgl,filebase+"_genelist")

# This is the workflow that creates the average (or "metagene") plots.
# alignpos is a variable: 0 means align 5'end, 1 means start codon, 2 means stop codon, 3 means 3' end of mRNA as defined with bp5, bp3.
# regionlength5 and 3 are the distances up and down of the alignment position that are considered in the average. Note that UTR inputs need to be -1 or the bp distances (+UTR lengths) big enough to accomodate the regionlengths.

# This is for creating an average of many genes.
# Inputs are:
# regionlength - distance over which to average genes (minimal length of gene).
# countlist - list of all genes to be averaged. 
# equalweight - set to 1 if all genes are to be weighted equally in average, otherwise 0.
# One concern with this is still that genes with missing reads due to mapping to more than 1 location could influence this.

def totalavg_wf(regionlength5,regionlength3,filebase,counts_filestring,bp5,bp3,
				ignoreutr5,ignoreutr3,shift,filter,thresh,equalweight,
				GFFgen_filename,utrgfffilename,utr5gfffilename,alignpos,totalreads):
	#comments="This is a test run of totalavg_wf.\nIt makes an average of all genes.\n"
	#comments+="shift is amount reads are shifted for calculation of gene counts in the list.\n"
	#comments+="bp5 and bp3 are the extra distance on 5 and 3' ends for average. It is the entire UTR if no UTR info available for that gene. -1 for UTR means that all UTRs will be bp distance.\n"
	#comments+="thresh signifies minimal rpkm needed in coding region for gene to be in the average.\n"
	#comments+="alignpos =0 anchors average around the 5'UTR end and only includes 5'UTRs if UTR info given, otherwise illegal genes. If no UTR GFF, then it uses bp5. alignpos =3 is same for 3'UTR."
	#comments+="alignpos =1 anchors average around the start codon and only includes 5'UTRs if UTR info given, otherwise illegal genes. If no UTR GFF, then is uses bp5. alignpos =2 is same for stop codon."

	# Write output file of comments.
	#fc=open(filebase+"_"+str(alignpos)+"_output.txt","w")
	#fc.write(comments)
	#fc.write("\n")
	#fc.write("Totalavg_wf was called with parameters:\n")
	#fc.write("counts_filestring="+str(counts_filestring)+"\n")
	#fc.write("bp5,bp3="+str(bp5)+", "+str(bp3)+"\n")
	#fc.write("shift="+str(shift)+"\n")
	#fc.write("filter="+str(filter)+"\n")
	#fc.write("thresh="+str(thresh)+"\n")
	#fc.write("regionlength5="+str(regionlength5)+"\n")
	#fc.write("regionlength3="+str(regionlength3)+"\n")
	#fc.write("equalweight="+str(equalweight)+"\n")
	#fc.write("alignpos="+str(alignpos)+"\n")
	#fc.write("ignoreutr5="+str(ignoreutr5)+"\n")
	#fc.write("ignoreutr3="+str(ignoreutr3)+"\n")
	#fc.close()

	f=open(filebase+"_avg_"+str(alignpos),"wb")

	GFFgen=GFF.parse(GFFgen_filename)
	GFFlist=seqtools.makeGFFlist(GFFgen)

	counts0=builddense.readcountsf(counts_filestring+"_plus_")		
	counts1=builddense.readcountsf(counts_filestring+"_minus_")
	counts=[counts0,counts1]
	
	if ignoreutr3 != 1:
		utrgffgen = GFF.parse(utrgfffilename)
		utrtable = genometools.makeutrtable(utrgffgen)
	else:
		utrtable = {}
		
	if ignoreutr5 != 1:
		utrgffgen = GFF.parse(utr5gfffilename)
		utrtable2 = genometools.makeutrtable(utrgffgen)
	else:
		utrtable2 = {}
		
	gene=makeavggene(regionlength5,regionlength3,counts,GFFlist,utrtable,utrtable2,
					ignoreutr5,ignoreutr3,bp5,bp3,shift,filter,
					thresh,alignpos,equalweight,totalreads)
	
	#output dict of avgvalues to bokeh in avggene_wf
	values = {}
	for i in range(0,regionlength5+regionlength3):
		f.write(struct.pack("f",float(gene[i])))
		values[i] = float(gene[i])
		
	f.close()
	return values
	

	
# Function to write dict out to csv format.
# Modified on 5/7/12 to be backward compatible, but to also have the new feature of being able to 
# write out many dicts sequentially to the same csv.
def writedicttoexcel(genelists,filestring):
	import csv
	writer = csv.writer(open(filestring+".csv", "wb"),delimiter=',')
	if type(genelists)!=list:
		genelists=[genelists]		### This is also new code -- converting any input dict to a length 1 list.
	for genelist in genelists:		### This is the new code (just this line and changing genelist to genelists in func def.)
		
		if(genelist.has_key("headers")):		# This here presumably to get header as 1st line in csv.
			headerrecord=[]
			headerrecord.append("headers")
			for field in genelist["headers"]:
				headerrecord.append(field)
			writer.writerow(headerrecord)
		for gene in genelist.keys():
			generecord=[]
			generecord.append(gene)
			
			# New Feb 2013, check to see if we have a list or a single value.
			if type(genelist[gene])==list:
				for field in genelist[gene]:
					generecord.append(field)
			else:
				generecord.append(genelist[gene])
			if gene=="headers":				# Skip since we did this above.
				continue
			writer.writerow(generecord)
		
# Make a csv back into a dictionary. f is open file handle. Handles multiple hits for same gene by putting in a _# term for them.
def readindict(f):
	previousgene=""
	counter=1
	filegen=csv.reader(f,delimiter=',')
	output = {}
	for gene in filegen:
		if gene[0]==previousgene:
			modgenename=gene[0]+"_"+str(counter)
			counter+=1
		else:
			modgenename=gene[0]
			counter=1
		output[modgenename]=[]
		for column in gene[1:]:
			output[modgenename].append(column)
		previousgene=gene[0]
	f.close()	
	return output


# Convert makelist sequence output to a fasta file just writing out the incoming list of strings.
# genelist is the dictionary keyed to gene id, trans is set to 2 if the sequences should be translated. f_outfile is handle to writeable file.
# set trans to 0 if not translated and 1 for no trnaslation but conversion of * to X.
def makefastas(f_outfile,genelist,trans):
	for item in genelist.keys():
		if item=="headers":
			continue		
		fastaseq=SeqIO.SeqRecord(Seq(genelist[item][-1]))
		if str(fastaseq)=="":
			continue		
		if trans>=1:
			fastaseq=fastaseq.seq
		if trans==2:
			fastaseq=fastaseq.translate()
		if trans>=1:
			# Stop codons will be replaced with "X"
			for i in range(len(fastaseq)):
				if fastaseq[i] == '*':
					fastaseq=fastaseq[0:i]+Seq("X")+fastaseq[i+1:]
			fastaseq=SeqRecord(fastaseq)
		fastaseq.id=item
		f_outfile.write(fastaseq.format("fasta"))
	f_outfile.close()


# This workflow will be used to take in a csv that was modified in excel (to get the top readthrough candidates).
# It will write out the fastas of those stop codon contexts.
# outfilebase is where the fasta should go (path and filebase), infile is path and filename of input file.
# No official comment here since this is a simple function -- notes should be recorded though whenever used.
# Set prot to 1 for conversion of * to X. 0 for dna. 2 for translation.
def fastaout_wf(outfilebase,infile,prot):
	if prot>=1:
		f_outfile=open(outfilebase+"_prot.fasta","w")
	else:
		f_outfile=open(outfilebase+"_DNA.fasta","w")
	f=open(infile)
	genelist=readindict(f)
	makefastas(f_outfile,genelist,prot)
	f_outfile.close()
	f.close()
			
# This function takes in a csv that is generated by listavg and computes averages for particular columns dependent on other columns.
# It is intended to work with termination parameters.
# It then takes in 2 lists:
# list1 is a list of the columns containing the data values to be averaged.
# list2 is a list of the columns containing parameters to select which values are averaged in list1.
# The output is a a outfilebase where the new csv will be put. 
def termprops_wf(infilename,outfilename,list1,list2):
	f=open(infilename)
	indict=readindict(f)	
	
	# Eliminate neg values.
	negkeys=[]
	for col in list1:
		for keyval in indict:
			if indict[keyval][col][0]=='-':
				negkeys.append(keyval)
	if len(negkeys)>0:
		for keyval in negkeys:
			del indict[keyval]	
		print "Negative values detected and eliminated from dictionary."
		
	genelists=[]
	for item in list2:
		genelists.append({})
	
	headerlist=indict["headers"]
	# Make dicts.
	for gene in indict.keys():
		if gene=="headers":				
		# Take care of headers for each dict specified in list2.	
			list2pos=0	 # position we are at in list2. 
			for column in list2:
				genelists[list2pos]["headers"]=[]
				list1pos=0	 # which item we are at in list1.
				for col in list1:
					genelists[list2pos]["headers"].append(indict["headers"][col]+"_"+indict["headers"][column])
					list1pos+=1
				genelists[list2pos]["headers"].append("Number")		# This will be the count of each.
				list2pos+=1
			continue
			
			
		list2pos=0		 # position we are at in list2.
		for column in list2:		
			# Compute each average.
			
			# Check if we've used this key yet:
			# If not, create a list with a list of length of list1 with 0s.
			if not genelists[list2pos].has_key(indict[gene][column]):
				genelists[list2pos][indict[gene][column]]=[float(0) for x in range(len(list1)+1)]	# The + 1 here is so the last element can be a counter of the occurence of each key.
			
			list1pos=0	# which item we are at in list1.
			for col in list1:					
				genelists[list2pos][indict[gene][column]][list1pos]+=float(indict[gene][col])	# This will not hit last counter element.
				list1pos+=1	
			genelists[list2pos][indict[gene][column]][list1pos]+=1			#Now update counter element.
			list2pos+=1
			
	# Now compute averages:
	list2pos=0
	for column in list2:
		for entry in genelists[list2pos].keys():
			if entry=="headers":
				continue
			list1pos=0
			for col in list1:
				
				# Skip headers
				
				genelists[list2pos][entry][list1pos]/=genelists[list2pos][entry][len(list1)]
				list1pos+=1
		list2pos+=1
	writedicttoexcel(genelists,outfilename)			
	f.close()	
			
			

# Genelist creator and tools for exporting/importing. This is still going to be a complex, workhorse function as it was in the long past. But now it outputs sequences and uses givegene, so it's better. No rtrunc anymore either.
# This function is called by above workflow.
# Note that as of 2/9/13, it no longer outputs the dictionary of counts for each gene -- a new function will do that called makeavggene.
# In the future, it would make sense to make similar functions for uORFs perhaps.
# And the search functions that compare ratios at basepair level, will also do this sort of thing and output the same data.
# bp input is a parameter that extends the natural UTRs.
# shift is the distance we assume between the 5' end of the ribosome and the base it is at. 0 can be used for mRNA-seq on ORFs.
# filter is T or F, do you need extra filtering
# The output is a dictionary of gene info, summed counts, sequences, and other info. Dictionaries are keyed by gene id.
# thresh is the minimal number of normalized reads allowed on a gene -- in rpkm terms.
# The input counts is a list of 2 counts lists actually, strand 1 being the first, -1 the second.
# UPDATE May 2012: this function will now also output a second extra field (the 1st extra is readthrough). extra2 will be a termination pause score.
# It will also output some sequence info: penultimate aa, ultimate aa, 1st aa after, stop codon, ultimate bp, 1st bp after. Of course full sequence is still there too.
# May 8, 2012: also added tetra nucleotide motif and post-post bp.
# May 27, 2012: The old readthrough of just "bp" bases into 3' UTR is commented out and instead it takes in the utrtable for the 3'UTR and uses that.
# Note that putting in a utrignore of 1 will do it using bp inputs. This would be the way to do it for creating averages without 3'UTR as well since counts will be of 3'UTRs only otherwise.
# Make sure bp3_0 and bp5_0 are 0 if you want ONLY ORFs.
# bp has also now been changed to be 2 parameters for 5' and 3' ends. To get 3'UTR length, take length of output of sequence of 3'UTR.
# Bug found on 10/30/12. Genes without 3'UTR annotation weren't always being with no 3'UTR. Now fixed.
# Added on 10/30/12: Outputs for extrareads (Readthrough) for 2ndhalf of 3'UTR only (UTR length defined as annotated 3'UTR plus bp3). Also, peak position in 3'UTR and 3'UTR sequence.
# Added on 11/15/12 is 5'UTR sequence as extrasequence2. This comes in as utrtable2.  
# And now extra reads are reported as rpkm, that is not as fraction of main ORF.
# Note that 1 can be put in for ignoreutr5 or ignoreutr3 to not consider them in computing the "gene" and instead just use the bp3 and bp5 values.

def makegenelist(counts,GFFlist,utrtable,utrtable2,bp5_0,bp3_0,ignoreutr5,ignoreutr3,
				shift, filter, thresh, totalreads):
	genelist={}		# The list we will output.
	missedthresh=0
	illegalgenes=0
	genesinlist=0
	
	if bp5_0<0 or bp3_0<0:
		print "Illegal negative values for bp."
		exit()
	
	# Put headers on genelist
	genelist["headers"]=["alias","feat_num","rpkm", "rpc"]
	for chrom in GFFlist:
		feat_num=0
		for feature in GFFlist[chrom].features:
			if ignoreutr3!=1:
				if utrtable.has_key(feature.id):
					bp3=utrtable[feature.id][1]-utrtable[feature.id][0]+bp3_0
					noutr3=False
				else:
					bp3=bp3_0	
					noutr3=True
			else:
				bp3=bp3_0
				noutr3=False
			
			if ignoreutr5!=1:
				if utrtable2.has_key(feature.id):
					bp5=utrtable2[feature.id][1]-utrtable2[feature.id][0]+bp5_0
					noutr5=False
				else:
					bp5=bp5_0	
					noutr5=True
			else:
				bp5=bp5_0
				noutr5=False
					
			# Import sequence and counts. This is using the new (as of Feb 2013) givegene which takes shift as an input parameter.
			gg=seqtools.givegene(chrom,feat_num,[GFFlist,utrtable2,utrtable],counts,[bp5,bp3,shift],2)
			genesequence=gg[1]
			genecounts=gg[0]

			# Get rid of dubious genes, nongenes, genes with overlap of others.
			if (genesequence==-1 or genecounts ==-1 or genesequence==-2 or genecounts==-2):
				feat_num+=1
				illegalgenes+=1
				continue
				
			# Define ORF
			start=bp5
			end=len(genecounts)-bp3
			if start==end:
				print "Error, gene length is 0 for gene "+feature.id
				exit()
				
			# Compute rpkm for each.
			totalgenereads=float(1000)*sum(genecounts[start:end])/len(genecounts[start:end])

			# compute reads per codon, added 7/21/2014 ab
			# threshold on this now !! not rpkm
			rpm_per_codon = 3 * sum(genecounts[start:end]) / len(genecounts[start:end])
			readfactor = totalreads / float(1000000)
			readspercodon = readfactor * rpm_per_codon
			if readspercodon < thresh:
				feat_num += 1
				missedthresh += 1
				continue
					
			# Filter sequences for those of interest - Here is where the action is. i.e. specific sequences, or go term.
			if filter:
				has_seq = filtermodule(feature, GFFlist, chrom)
				if not has_seq:
					feat_num+=1
					continue		
				
			if "Alias" in feature.qualifiers:
				alias = feature.qualifiers["Alias"][0]
			elif "Name" in feature.qualifiers:				# For coli - added 5/9/12.
				alias=feature.qualifiers["Name"][0]
			else:
				alias = "NA"
						
			if "Note" in feature.qualifiers:
				note = feature.qualifiers["Note"][0]
			else:
				note = "NA"		
		
			# Output chrom, featurenum, counts, alias.

			genelist[feature.id]=[]
			genelist[feature.id].append(alias)
			genelist[feature.id].append(feat_num)	
			genelist[feature.id].append(totalgenereads)	
			genelist[feature.id].append(readspercodon)

			feat_num+=1	
			genesinlist+=1

	print "Genes below threshold = " + str(missedthresh)
	print "Genes dropped by givegene (overlap, undesirable features, etc.) = " + str(illegalgenes)
	print "Genes in list = " + str(genesinlist)
	
	return genelist
	

# These Properties of gene list that always will be true:
# Keyed on gene id
# First column alias, second chrom id, third featurenum, last is sequence and penultimate is note.
# In between are reads, with the most important or interesting coming first. Other info comes after that.
# There is also a header key called "header".			
			

# This function simply averages the counts for every gene. It borrows much from makegenelist.
# For now, I don't think it makes sense to make a common GFF loop program, as the needs are a bit different here. Eventually though that may be the way to go.
# alignpos is a variable: 0 means align 5'end, 1 means start codon, 2 means 1 bp past stop codon, 3 means 1 bp past the 3' end of mRNA as defined with BOTH the UTR annotation and bp5_0, bp3_0.
# Note that alignpos = 0 or 3 means that the ORF will not be included in average! And in this case, genes without UTR annotation are excluded unless ignoreutr inputs are 1.
# If ignoreutr is 0, then any alignment = 1 or 2 will require that the UTR and ORF are long enough to accomodate the regionlengths. 
# bp5_0 and bp3_0 are extra distances that will be added to the annotated 5'and 3' UTR lengths.
# regionlength5 and 3 are the distances up and down of the alignment position that are considered in the average.
# so average plot is length regionlength5+regionlength3, with alignent position the 1st in regionlength3.
# This function is called by above workflow.
# Note that alignpos0,3 will maintain the reading frame by adding +-1 nt to bp3_0,bp5_0 to keep it in frame with ORF.
# One concern with this is still that genes with missing reads due to mapping to more than 1 location could influence this.
def makeavggene(regionlength5,regionlength3,counts,GFFlist,utrtable,utrtable2,
				ignoreutr5,ignoreutr3,bp5_0,bp3_0,shift,filter,
				thresh,alignpos,equalweight,totalreads):
	missedthresh=0
	illegalgenes=0
	genesinlist=0
	tooshortlist=0
	averagegene=[0 for x in range(0,(regionlength5+regionlength3))]

	if bp5_0<0 or bp3_0<0:
		print "Illegal negative values for bp."
		exit()
	
	for chrom in GFFlist:
		feat_num=-1
		for feature in GFFlist[chrom].features:
			feat_num+=1
			if ignoreutr3!=1:
				if utrtable.has_key(feature.id):
					bp3=utrtable[feature.id][1]-utrtable[feature.id][0]+bp3_0
					noutr3=False
				else:
					bp3=bp3_0	# This line doesn't do anything because noutr=True genes are dumped.
					noutr3=True
			else:
				bp3=bp3_0
				noutr3=False
			
			if ignoreutr5!=1:
				if utrtable2.has_key(feature.id):
					bp5=utrtable2[feature.id][1]-utrtable2[feature.id][0]+bp5_0
					noutr5=False
				else:
					bp5=bp5_0	# This line doesn't do anything because noutr=True genes are dumped.
					noutr5=True
			else:
				bp5=bp5_0	
				noutr5=False
			
			
			# Get proper extensions on gene. Check for cases here where UTR will be too short.
			if alignpos==0:		# Align on 5' ends.
				if noutr5==True:
					illegalgenes+=1
					continue
				elif regionlength3>bp5:
					tooshortlist+=1
					continue
				else:	
					bp5=regionlength5+bp5
					frshift=bp5%3
					if frshift==2:
						frshift=-1
					bp5-=frshift


			elif alignpos==1:	# Align at start codon.
				if noutr5==True:
					illegalgenes+=1
					continue
				elif regionlength5>bp5:
					tooshortlist+=1
					continue
				else:
					bp5=regionlength5
					
			elif alignpos==2:	#align at stop codon.
				if noutr3==True:
					illegalgenes+=1
					continue
				elif regionlength3>bp3:
					tooshortlist+=1
					continue
				else:
					bp3=regionlength3
					
			elif alignpos==3:	# Align at 3' end.
				if noutr3==True:
					illegalgenes+=1
					continue
				elif regionlength5>bp3:
					tooshortlist+=1
					continue
				else:
					bp3=regionlength3+bp3
					frshift=bp3%3
					if frshift==2:
						frshift=-1
					bp3-=frshift
			
			# Get counts
			gg=seqtools.givegene(chrom,feat_num,[GFFlist,utrtable2,utrtable],counts,[bp5,bp3,shift],2)
			genecounts=gg[0]
			genesequence=gg[1]
			# Get rid of dubious genes, nongenes, genes with overlap of others.
			if (genesequence==-1 or genecounts ==-1 or genesequence==-2 or genecounts==-2):
				illegalgenes+=1
				continue

			# Define ORF
			start=bp5
			end=len(genecounts)-bp3
			if start==end:
				print "Error, gene length is 0 for gene "+feature.id
				exit()
						
			# Filter sequences for those of interest - Here is where the action is. i.e. specific sequences, or go term.
			if filter:
				has_seq = filtermodule(feature, GFFlist, chrom)
				if not has_seq:
					continue	
			
# 			## Threshold on total gene counts, rpkm 
# 			totalgenereads=float(1000)*sum(genecounts[start:end])/len(genecounts[start:end])
# 			if totalgenereads<thresh:
# 				missedthresh+=1
# 				continue
	
			# Threshold on reads per codon, added 7/18/2014 ab
			rpm_per_codon = 3 * sum(genecounts[start:end]) / len(genecounts[start:end])
			readfactor = totalreads / float(1000000)
			reads_per_codon = readfactor * rpm_per_codon
			if reads_per_codon < thresh:
				missedthresh += 1
				continue
	
			#Cut down to the size we want and make sure we don't run off the rails.
			if alignpos==0:
				# All checks for length already done above.
				genesinlist+=1
				countlist=genecounts[0:(regionlength5+regionlength3)]
				if equalweight==1:
					totcounts=sum(countlist)
				else:
					totcounts=1
				# Add to growing average.
				if totcounts!=0:	
					for i in range(len(countlist)): 
						averagegene[i]+=countlist[i]/totcounts
			
			elif alignpos==1:
				# Check for enough length on 3' end (5' end already checked above):
				if len(genecounts)<(regionlength5+regionlength3):
					tooshortlist+=1
					continue
				else:
					genesinlist+=1
					countlist=genecounts[0:(regionlength5+regionlength3)]
					if equalweight==1:
						totcounts=sum(countlist)
					else:
						totcounts=1
					# Add to growing average.
					if totcounts!=0:	
						for i in range(len(countlist)): 
							averagegene[i]+=countlist[i]/totcounts
			
			
			elif alignpos==2:
				# Check for enough length on 5' end (3' end already checked above):
				if len(genecounts)<(regionlength5+regionlength3):
					tooshortlist+=1
					continue
				else:
					genesinlist+=1
					countlist=genecounts[-(regionlength5+regionlength3):]
					if equalweight==1:
						totcounts=sum(countlist)
					else:
						totcounts=1
					# Add to growing average.
					if totcounts!=0:	
						for i in range(len(countlist)): 
							averagegene[i]+=countlist[i]/totcounts
							
							
			elif alignpos==3:
				# All checks for length already done above.
				genesinlist+=1
				countlist=genecounts[-(regionlength5+regionlength3):]
				if equalweight==1:
					totcounts=sum(countlist)
				else:
					totcounts=1
				# Add to growing average.
				if totcounts!=0:	
					for i in range(len(countlist)): 
						averagegene[i]+=countlist[i]/totcounts
			
	if equalweight==0:
		m=0
		for m in range(len(averagegene)):
			if genesinlist!=0:
				averagegene[m]/=genesinlist
			else:
				print "Error, no genes to average at start."		
					
	print "Genes under threshold: "+str(missedthresh)
	print "Genes removed for overlap or no UTR or dubious or non-major chrom or other: "+str(illegalgenes)
	print "Genes removed because feature(s) too short: "+str(tooshortlist)
	print "Genes in average: "+str(genesinlist)
	
	return averagegene


## tests to see if the feature contains some sequence element, returns True or False
def filtermodule(feature, GFFlist, chrom):
	motif1 = 'PP'
	geneseq = feature.extract(GFFlist[chrom]).seq
	geneseq = geneseq.translate()
	if motif1 in geneseq:
		test = True
	else:
		test = False
	return test


# 	# A module to filter by GO term.
# 	# Currently set to select for membrane proteins.
# def GOfilt(feature):
# 
# 	for term in feature.qualifiers["Ontology_term"]:
# 		if term == "GO:0016020": # Membrane proteins.
# 			return 1
# 	return 0

	
# A Function to make a masterdict
# Change on 2/16/13 - do not put in sequence for full chromosomes.
def makemasterdict(GFFgen):
	masterdict={}
	masterdict["headers"]=["alias","chromosome","featurenum","Note","Sequence"]
	for chr in GFFgen:
		feat_num=0
		for feature in chr.features:
			if "Alias" in feature.qualifiers:
				alias = feature.qualifiers["Alias"][0]
			elif "Name" in feature.qualifiers:	
				alias=feature.qualifiers["Name"][0]			
			else:
				alias = "NA"
			if "Note" in feature.qualifiers:
				note = feature.qualifiers["Note"][0]
			else:
				note = "NA"
			start=feature.location.start.position
			end=feature.location.end.position
			sequence=chr[start:end].seq
			if feature.strand==-1:
				sequence=sequence.reverse_complement()
			seq=str(sequence)
			if feat_num==0:
				seq=""		# So we don't have whole chromosome here.							
			masterdict[feature.id]=[alias,chr.id,feat_num,note,seq]
			feat_num+=1	
	return masterdict
	
# A function to add a dictionary to a master dictionary for all genes.
# appendcollist is a list of the columns in newdict to append in every time.
# namelist is the names of those columns
# This works with any dictionary where the dictionary (newdict) has the keys of feature id.
def combinetomaster(masterdict,newdict,GFFgen,appendcollist,namelist):
	if len(appendcollist)!=len(namelist):
		print "error in masterlist combiner."
	for colname in namelist:
		masterdict["headers"].append(colname)
	for chr in GFFgen:
		for feature in chr.features:
			for colnum in appendcollist:
				if newdict.has_key(feature.id):
					masterdict[feature.id].append(newdict[feature.id][colnum])
				else:
					masterdict[feature.id].append(-10)
				
	return masterdict
