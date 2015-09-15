import ConfigParser
import sys
import subprocess
import time
from collections import defaultdict


def alpha(int_list):
	ans = ['0']*len(int_list)
	mapper="0ABCDEFGHIJKLMNPQRSTUVWXYZ"
	for pos in range(len(int_list)):
		ans[pos] = mapper[int_list[pos]]
	return ans

def prepMember(member, configpath):
	"""
	Assuming configpath is a config file and member is either Child,
	Mother of Father. f1 and f2 are the paths to the files where the
	two sequences for this member are going to be scored.
	"""
	config = ConfigParser.RawConfigParser()
	config.read(configpath)
	ref = config.get("Common", "reference")
	win = config.get("Common", "window")
	vcf = config.get(member, "vcf_file")
	f1 = config.get(member, "fasta1")
	f2 = config.get(member, "fasta2")
	return [ref, win, vcf, f1, f2]

def doesLineCount(line):
	"""
	Assuming that line is a line from a VCF(v4.1) file, this function
	returns a True/False value depending on whether or not we should use
	this line of the VCF file. The function returns True if and only if
	the variant is present and is an insertion, deletion or replacement.
	"""
	if (
		line.find("[")==(-1) and
		line.find("]")==(-1) and
		line.find("<CGA_CNVWIN>")==(-1) and
		line.find("IMPRECISE")==(-1) and
		line.find("<CGA_NOCALL>")==(-1) and
		not (line.split())[4] == "."
	):
		return True
	else:
		return False

def VCFtoFASTA(member):
	"""
	Assuming that member gives us the information about a child, a mother
	or a father of the trio. We first build the reference window from the FASTA
	file, then apply the variants from the VCF file to form two FASTA sequences
	-- the first one takes all non-overlapping variants from left to right, and
	the second one takes the ones that can't be put in the first, as well as the
	homozygous variants. This process is described in [Makinen and Valenzuela, 2015].

	After building the two sequences, we align them by putting gaps where insertions
	and deletions have occured. This is not optimal alignment, but it works for us.
	Then, we write these aligned sequences to the FASTA files.
	"""
	ref = open(member[0], "r")
	window = open(member[1], "r")
	vcffile = open(member[2], "r")
	fasta1 = open(member[3], "w")
	fasta2 = open(member[4], "w")
	
	chrom = ">"+window.readline() #in the FASTA genome file, every chromosome is a header and starts with ">"
	start = int(window.readline()) #this is the start of the window (included in the reference window)
	finish = int(window.readline())-1 #this is the end of the window (not included in the reference window)

	refwin = [] #a list that corresponds to the reference window between the start and the end


	#first, we get the FASTA sequence between positions start and finish
	chrnum = ref.readline()

	while not chrnum == chrom: #make sure we're at the right chromosome
		chrnum = ref.readline()

	posskipped = 0

	for line in ref:
		line_length = len(line.strip())
		posskipped+=line_length #this is how long we assume that a FASTA line is
		if posskipped<start:
			continue
		refwin+=list(line.strip().upper())
		if posskipped>finish:
			break
	refwin = refwin[(start-1)%line_length:] #trim chars before the start position
	refwin = refwin[:-(line_length-finish%line_length)] #trim chars after the finish; turns out that if finish%line_length==0, then refwin has +line_length chars

	sequence1 = list(refwin) 
	sequence2 = list(refwin) 
	used1 = [0]*len(refwin) #indicates if an index is used to check for overlapping variants. It uses variant_id to asociate the position to a specific variant.
	used2 = [0]*len(refwin) #see above
	ind1 = [0]*len(refwin) #indicates where are insertions/deletions
	ind2 = [0]*len(refwin)
	var_map1 = [0]*len(refwin)
	var_map2 = [0]*len(refwin)

	homozygous = False #a homozygous variant is applied to both sequences

	#second, we insert the changes from the VCF into the sequences
	variant_counter = 0
	hetero_counter = 0

	for line in vcffile:
		if line[0]=='#': #header line
			continue
		if not doesLineCount(line):
			continue
		info = line.split()
		pos = int(info[1])
		if not info[0] == chrom[1:-1]: #wrong chromosome
			continue
		
		record_id = int(info[2].strip())
		ref_seq = list(info[3].strip())
		alt_seq = list(info[4].strip())
		if pos<start: #haven't reached the window yet
			continue
		if pos>finish: #went past the window
			break
		if pos + len(ref_seq) > finish: #variant will go out of the window
			break
		if ',' in alt_seq:
			continue #we don't want to see the commas for now TODO: divide into two lines, check VCF to make sure these aren't phased
		#we are surely in the window
		#now we check whether the variant is homozygous or heterozygous
		
		#print "processing the following line:"+str(record_id)
		#print line
		variant_counter = variant_counter + 1;
		variant_id = record_id 
		homozygous = True if info[9][0]=='1' and info[9][2]=='1' else False
		if homozygous == False:
			hetero_counter += 1
		winpos = pos-start
		assert(ref_seq == refwin[winpos:winpos+len(ref_seq)])
		offset_1 = sum(ind1[:winpos]) #what is the difference in positions between ref and sequence1
		offset_2 = sum(ind2[:winpos]) #what is the difference in positions between ref and sequence2
		uscore1 = sum(used1[winpos:winpos+len(ref_seq)]) #is there anything used in the positions of the variant
		uscore2 = sum(used2[winpos:winpos+len(ref_seq)])
		if uscore1 == 0: #we can add this variant to the 1st FASTA file	
			sequence1[winpos+offset_1:winpos+len(ref_seq)+offset_1] = list(alt_seq)
			var_map1[winpos+offset_1:winpos+len(ref_seq)+offset_1] = [variant_id]*len(alt_seq)
			ind1[winpos] = len(alt_seq)-len(ref_seq)
			used1[winpos:winpos+len(ref_seq)] = [variant_id]*len(ref_seq) # val uses ref_seq. alt_seq keeps used1 aligned to seequence1
			assert(variant_id > 0)
			if not homozygous: #we wouldn't want to add to the second sequence
				continue
		if uscore2 == 0:
			sequence2[winpos+offset_2:winpos+len(ref_seq)+offset_2] = list(alt_seq)
			var_map2[winpos+offset_2:winpos+len(ref_seq)+offset_2] = [variant_id]*len(alt_seq)
			ind2[winpos] = len(alt_seq)-len(ref_seq)
			used2[winpos:winpos+len(ref_seq)] = [variant_id]*len(ref_seq) # same as above
			assert(variant_id > 0)

	#print "We processed :"+str(variant_counter)+" variants"
	#print "We processed :"+str(hetero_counter)+" heterozygous variants"
	#print "********"
	#print "".join(map(str,used1))
	#print "".join(map(str,used2))
	#print "********"
	#print "".join(sequence1)
	#print "".join(sequence2)
	#print "********"
	#print "".join(map(str,ind1))
	#print "".join(map(str,ind2))
	#print "********"
	gappedsequence1 = list(sequence1) #in these we put the aligned sequences
	gappedsequence2 = list(sequence2)
	happening = 0
	gaps1 = 0
	gaps2 = 0
	#now, we build the alignment of the sequences (in fact, only semi-align: put gaps where there are INDELS)
	#put gaps ('-'s) with respect to the first sequence
	for nucind in xrange(len(ind1)):
		gaplen = ind1[nucind]*((-1) if ind1[nucind]<0 else 1)
		if ind1[nucind]>0: #there was an insertion
			#add gaps to sequence2 and make sure ind is OK too
			gappedsequence2[nucind+gaps2:nucind+gaps2+1] = gappedsequence2[nucind+gaps2:nucind+gaps2+1] + ['-']*gaplen
			var_map2[nucind+gaps2:nucind+gaps2+1] = var_map2[nucind+gaps2:nucind+gaps2+1] + [0]*gaplen
			
			gaps1+=gaplen
			gaps2+=gaplen
			happening+=1
		if ind1[nucind]<0: #there was a deletion
			#add gaps to sequence1 and make sure ind is OK too
			gappedsequence1[nucind+gaps1:nucind+gaps1+1] = gappedsequence1[nucind+gaps1:nucind+gaps1+1] + ['-']*gaplen
			var_id = var_map1[nucind+gaps1]
			assert(var_id != 0)
			var_map1[nucind+gaps1:nucind+gaps1+1] = var_map1[nucind+gaps1:nucind+gaps1+1] + [var_id]*gaplen  #TODO
			happening+=1
	# put gaps with respect to the second sequence
	# NEEDS TO BE DONE IN THE SAME CYCLE TO KEEP THE GAPS CORRECT
		gaplen = ind2[nucind]*((-1) if ind2[nucind]<0 else 1)
		if ind2[nucind]>0: #there was an insertion
			#add gaps to sequence1 and make sure ind is OK too
			gappedsequence1[nucind+gaps1:nucind+gaps1+1] = gappedsequence1[nucind+gaps1:nucind+gaps1+1] + ['-']*gaplen
			var_map1[nucind+gaps1:nucind+gaps1+1] = var_map1[nucind+gaps1:nucind+gaps1+1] + [0]*gaplen
			
			gaps1+=gaplen
			gaps2+=gaplen
			happening+=1
		if ind2[nucind]<0: #there was a deletion
			#add gaps to sequence2 and make sure ind is OK too
			gappedsequence2[nucind+gaps2:nucind+gaps2+1] = gappedsequence2[nucind+gaps2:nucind+gaps2+1] + ['-']*gaplen
			var_map2[nucind+gaps2:nucind+gaps2+1] = var_map2[nucind+gaps2:nucind+gaps2+1] + [0]*gaplen
			happening+=1

	#lastly, we format the aligned sequences as FASTA (60 characters per line)
	fasta1.write(">hg19|chromosome "+chrom[1:-1]+"|start pos "+str(start)+"|end pos "+str(finish)+"|variants from "+member[2]+"|first sequence\n")
	fasta2.write(">hg19|chromosome "+chrom[1:-1]+"|start pos "+str(start)+"|end pos "+str(finish)+"|variants from "+member[2]+"|second sequence\n")
	fasta1.write("\n".join("".join(gappedsequence1[i:i+60]) for i in xrange(0, len(gappedsequence1), 60))+"\n")
	fasta2.write("\n".join("".join(gappedsequence2[i:i+60]) for i in xrange(0, len(gappedsequence2), 60))+"\n")
	#print "".join(gappedsequence1)
	#print "".join(gappedsequence2)
	#print "______"
	#print "".join(map(str,alpha(var_map1)))
	#print "".join(map(str,alpha(var_map2)))
	#print "______"

	vcffile.close()
	ref.close()
	window.close()

	fasta1.close()
	fasta2.close()

	#here we return the two ind lists
	#return [ind1, ind2]
	return [var_map1, var_map2, hetero_counter]

def callSimilarityPhaser(mother, father, child):
	mother_fasta1 = mother[3]
	mother_fasta2 = mother[4]

	father_fasta1 = father[3]
	father_fasta2 = father[4]

	child_fasta1 = child[3]
	child_fasta2 = child[4]

	print mother_fasta1, mother_fasta2, child_fasta1, child_fasta2
	#subprocess.check_call("make -C phasing_family/src", shell=True)
	quad_pass = 0;
	subprocess.check_call("phasing_family/src/mfc_similarity_phaser {0} {1} {2} {3} {4} {5} {6}".format(mother_fasta1, mother_fasta2, father_fasta1, father_fasta2, child_fasta1, child_fasta2, quad_pass), shell=True)



def phasedStringToVCF(child):
	var_map1 = child[5]
	var_map2 = child[6]
	hetero_vars_applied = child[7]
	stringfile = open("phase_string.txt", "r")
	phaseString = list(stringfile.read().strip())
	counts_dict = defaultdict(int)
	for pos in range(len(phaseString)):
		if (phaseString[pos] != '?'):
			switch = 0
			if phaseString[pos] == '1': # From father , 1|0
				switch = 1
			elif phaseString[pos] == '0': # From mother, 0|1
				switch = -1
			else:
				raise ValueError("UNKNOWN CHAR IN PHASE STRIGN, ABORTING")
			var_id_1 = var_map1[pos]
			if (var_id_1 != 0):
				# either adds or decreases one.
				counts_dict[var_id_1] += switch
			
			var_id_2 = var_map2[pos]
			if (var_id_2 != 0):
				counts_dict[var_id_1] -= switch
	print "Dict contains: " + str(len(counts_dict)) + " entries"
	vcffile = open(child[2], "r")
	new_vcffile = open(child[2]+".new.vcf", "w")
	
	right_count = 0
	wrong_count = 0
	no_gt_count = 0
	for line in vcffile:
		if line[0]=='#': #header line
			new_vcffile.write(line)
			continue
		info = line.split()
		record_id = int(info[2].strip())
		if record_id not in counts_dict:
			new_vcffile.write(line)
			continue	
		# record_id in counts_dict:
		count = counts_dict[record_id]
		if count > 0:
			phase = "1|0"
			no_phase = "0|1"
		elif count < 0:
			phase = "0|1"
			no_phase = "1|0"
		else:
			assert(count == 0)
			new_vcffile.write(line)
			continue	
		info = line.split('\t')
		real_phase = info[9][0:3]
		if real_phase == phase:
			right_count += 1
		elif real_phase == no_phase:
			wrong_count += 1
		else:
			no_gt_count += 1
		info[9] = phase
		new_line = '\t'.join(info)
		new_line += "\n"
		new_vcffile.write(new_line)
	print "-------------------------------------------" 
	print "Correct Phases:" +str(right_count)
	print "Wrong Phases:" +str(wrong_count)
	print "Phases with no ground truth info:" + str(no_gt_count) 
	print "-------------------------------------------" 
	print "Hetero Vars Applied:" +str(hetero_vars_applied)
	print "-------------------------------------------" 
def getFamilyFASTA():	
	try:	
		configuration = sys.argv[1]
		mother = prepMember("Mother", configuration)
		father = prepMember("Father", configuration)
		child = prepMember("Child", configuration)
		try: #attach the ind lists for the two FASTA files to the individual
			mother = mother +VCFtoFASTA(mother)
			father = father + VCFtoFASTA(father)
			child = child + VCFtoFASTA(child)
		except IOError:
			print "Please make sure the path to the resulting FASTA file is valid."
	except IndexError:
		print "Please select a configuration file."
	except:
		print "The configuration file you selected is invalid."
		raise
	return (mother, father, child)

(mother, father, child) = getFamilyFASTA()
callSimilarityPhaser(mother, father, child)
phasedStringToVCF(child)


