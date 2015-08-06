import ConfigParser
import sys
import subprocess

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
	finish = int(window.readline()) #this is the end of the window (included in the reference window)

	refwin = [] #a list that corresponds to the reference window between the start and the end


	#first, we get the FASTA sequence between positions start and finish
	chrnum = ref.readline()

	while not chrnum == chrom: #make sure we're at the right chromosome
		chrnum = ref.readline()

	posskipped = 0

	for line in ref:
		posskipped+=60 #this is how long we assume that a FASTA line is
		if posskipped<start:
			continue
		refwin+=list(line.strip().upper())
		if posskipped>finish:
			break
	refwin = refwin[(start-1)%60:] #trim chars before the start position
	refwin = refwin[:-(60-finish%60)] #trim chars after the finish; turns out that if finish%60==0, then refwin has +60 chars

	sequence1 = list(refwin) 
	sequence2 = list(refwin) 
	used1 = [0]*len(refwin) #indicates if an index is used to check for overlapping variants
	used2 = [0]*len(refwin) #see above
	ind1 = [0]*len(refwin) #indicates where are insertions/deletions
	ind2 = [0]*len(refwin)

	homozygous = False #a homozygous variant is applied to both sequences

	#second, we insert the changes from the VCF into the sequences
	for line in vcffile:
		if line[0]=='#': #header line
			continue
		if not doesLineCount(line):
			continue
		if not line[0] == chrom: #wrong chromosome
			continue
		info = line.split()
		pos = int(info[1])
		refpos = list(info[3].strip())
		altpos = list(info[4].strip())
		if pos<start: #haven't reached the window yet
			continue
		if pos>finish: #went past the window
			break
		#we are surely in the window
		#now we check whether the variant is homozygous or heterozygous
		homozygous = True if info[9][0]=='1' and info[9][2]=='1' else False
		winpos = pos-start
		si1 = sum(ind1[:winpos]) #what is the difference in positions between ref and sequence1
		si2 = sum(ind2[:winpos]) #what is the difference in positions between ref and sequence2
		uscore1 = sum(used1[winpos+si1:winpos+len(refpos)+si1]) #is there anything used in the positions of the variant
		uscore2 = sum(used2[winpos+si2:winpos+len(refpos)+si2])
		if uscore1 == 0: #we can add this variant to the 1st FASTA file
			sequence1[winpos+si1:winpos+len(refpos)+si1] = list(altpos)
			ind1[winpos] = len(altpos)-len(refpos)
			used1[winpos+si1:winpos+len(refpos)+si1] = [1]*len(refpos)
			if not homozygous: #we wouldn't want to add to the second sequence
				continue
		if uscore2 == 0:
			sequence2[winpos+si2:winpos+len(refpos)+si2] = list(altpos)
			ind2[winpos] = len(altpos)-len(refpos)
			used2[winpos+si2:winpos+len(refpos)+si2] = [1]*len(refpos)

	gappedsequence1 = list(sequence1) #in these we put the aligned sequences
	gappedsequence2 = list(sequence2)

	gaps1 = 0
	gaps2 = 0
	#now, we build the alignment of the sequences (in fact, only semi-align: put gaps where there are INDELS)
	#put gaps ('-'s) with respect to the first sequence
	for nucind in xrange(len(ind1)):
		gaplen = ind1[nucind]*((-1) if ind1[nucind]<0 else 1)
		if ind1[nucind]>0: #there was an insertion
			#add gaps to sequence2 and make sure ind is OK too
			gappedsequence2[nucind+gaps2:nucind+gaps2+1] = gappedsequence2[nucind+gaps2:nucind+gaps2+1] + ['-']*gaplen
			gaps2+=gaplen
		if ind1[nucind]<0: #there was a deletion
			#add gaps to sequence1 and make sure ind is OK too
			gappedsequence1[nucind+gaps1:nucind+gaps1+1] = gappedsequence1[nucind+gaps1:nucind+gaps1+1] + ['-']*gaplen
			gaps1+=gaplen	
	gaps1 = 0
	gaps2 = 0
	#put gaps with respect to the second sequence
	for nucind in xrange(len(ind2)):
		gaplen = ind2[nucind]*((-1) if ind2[nucind]<0 else 1)
		if ind2[nucind]>0: #there was an insertion
			#add gaps to sequence1 and make sure ind is OK too
			gappedsequence1[nucind+gaps1:nucind+gaps1+1] = gappedsequence1[nucind+gaps1:nucind+gaps1+1] + ['-']*gaplen
			gaps1+=gaplen
		if ind2[nucind]<0: #there was a deletion
			#add gaps to sequence2 and make sure ind is OK too
			gappedsequence2[nucind+gaps2:nucind+gaps2+1] = gappedsequence2[nucind+gaps2:nucind+gaps2+1] + ['-']*gaplen
			gaps2+=gaplen	
	#watch out: it might turn out that this alignment is not enough
	
	#lastly, we format the aligned sequences as FASTA (60 characters per line)
	fasta1.write(">hg19|chromosome "+chrom[1:-1]+"|start pos "+str(start)+"|end pos "+str(finish)+"|variants from "+member[2]+"|first sequence\n")
	fasta2.write(">hg19|chromosome "+chrom[1:-1]+"|start pos "+str(start)+"|end pos "+str(finish)+"|variants from "+member[2]+"|second sequence\n")
	fasta1.write("\n".join("".join(gappedsequence1[i:i+60]) for i in xrange(0, len(gappedsequence1), 60))+"\n")
	fasta2.write("\n".join("".join(gappedsequence2[i:i+60]) for i in xrange(0, len(gappedsequence2), 60))+"\n")
	
	vcffile.close()
	ref.close()
	window.close()
	
	fasta1.close()
	fasta2.close()
	
def callSimilarityPhaser(mother, father, child):
	mother_fasta1 = mother[3]
	mother_fasta2 = mother[4]
	
	father_fasta1 = father[3]
	father_fasta2 = father[4]
	
	child_fasta1 = child[3]
	child_fasta2 = child[4]
	
	print mother_fasta1, mother_fasta2, child_fasta1, child_fasta2
	subprocess.check_call("make -C phasing_family/src", shell=True)
	subprocess.check_call("phasing_family/src/mfc_similarity_phaser {0} {1} {2} {3} {4} {5} ".format(mother_fasta1, mother_fasta2, father_fasta1, father_fasta2, child_fasta1, child_fasta2), shell=True)
	
def vcfToFasta():	
	try:	
		configuration = sys.argv[1]
		mother = prepMember("Mother", configuration)
		father = prepMember("Father", configuration)
		child = prepMember("Child", configuration)
		try:
			VCFtoFASTA(child)
			VCFtoFASTA(mother)
			VCFtoFASTA(father)
		except IOError:
			print "Please make sure the path to the resulting FASTA file is valid."
	except IndexError:
		print "Please select a configuration file."
	except:
		print "The configuration file you selected is invalid."
	
	return (mother, father, child)

(mother, father, child) = vcfToFasta()
callSimilarityPhaser(mother, father, child)
#phasedStringToVcf()


