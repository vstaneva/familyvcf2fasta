import ConfigParser
import sys

def prepMember(member, configpath):
	config = ConfigParser.RawConfigParser()
	config.read(configpath)
	ref = config.get("Common", "reference")
	win = config.get("Common", "window")
	vcf = config.get(member, "vcf_file")
	f1 = config.get(member, "fasta1")
	f2 = config.get(member, "fasta2")
	return [ref, win, vcf, f1, f2]

def VCFtoFASTA(member):
	ref = open(member[0], "r")
	window = open(member[1], "r")
	vcffile = open(member[2], "r")
	fasta1 = open(member[3], "w")
	fasta2 = open(member[4], "w")

	start = int(window.readline())
	finish = int(window.readline())

	refwin = []


	#first, we get the FASTA sequence between positions start and finish
	ref.readline()

	posskipped = 0

	for line in ref:
		posskipped+=60 #this is how long a FASTA line is
		if posskipped<start:
			continue
		refwin+=list(line.strip().upper())
		if posskipped>finish:
			break
	refwin = refwin[(start-1)%60:] #trim chars before the start position
	refwin = refwin[:-(60-finish%60)] #trim chars after the finish; turns out that if finish%60==0, then refwin has +60 chars

	sequence1 = list(refwin)
	sequence2 = list(refwin)
	used1 = [0]*len(refwin)
	used2 = [0]*len(refwin)
	ind1 = [0]*len(refwin)
	ind2 = [0]*len(refwin)

	homozygous = False

	#second, we insert the changes from the VCF into the sequences
	for line in vcffile:
		if line[0]=='#':
			continue
		info = line.split()
		pos = int(info[1])
		refpos = list(info[3].strip())
		altpos = list(info[4].strip())
		if pos<start:
			continue
		if pos>finish:
			break
		#we are surely in that range
		#now we check whether the variant is homozygous or heterozygous
		homozygous = True if info[9][0]=='1' and info[9][2]=='1' else False
		winpos = pos-start
		si1 = sum(ind1[:winpos]) #what is the difference in positions between ref and sequence1
		si2 = sum(ind2[:winpos]) #what is the difference in positions between ref and sequence2
		uscore1 = sum(used1[winpos+si1:winpos+len(refpos)+si1])
		uscore2 = sum(used2[winpos+si2:winpos+len(refpos)+si2])
		if uscore1 == 0: #we can add this variant to the 1st FASTA file
			sequence1[winpos+si1:winpos+len(refpos)+si1] = list(altpos)
			ind1[winpos] = len(altpos)-len(refpos)
			used1[winpos+si1:winpos+len(refpos)+si1] = [1]*len(refpos)
			if not homozygous:
				continue
		if uscore2 == 0:
			sequence2[winpos+si2:winpos+len(refpos)+si2] = list(altpos)
			ind2[winpos] = len(altpos)-len(refpos)
			used2[winpos+si2:winpos+len(refpos)+si2] = [1]*len(refpos)

	gappedsequence1 = list(sequence1)
	gappedsequence2 = list(sequence2)

	gaps1 = 0
	gaps2 = 0

	#now, we build the alignment of the sequences

	for nucind in xrange(len(sequence1)):
		gaplen = ind1[nucind]*((-1) if ind1[nucind]<0 else 1)
		if ind1[nucind]>0: #there was an insertion
			#add gaps to sequence2 and make sure ind is OK too
			gappedsequence2[nucind+gaps2:nucind+gaps2+1] = gappedsequence2[nucind+gaps2:nucind+gaps2+1] + ['-']*gaplen
			gaps2+=gaplen
		if ind1[nucind]<0: #there was a deletion
			#add gaps to sequence1 and make sure ind is OK too
			gappedsequence1[nucind+gaps1:nucind+gaps1+1] = gappedsequence1[nucind+gaps1:nucind+gaps1+1] + ['-']*gaplen
			gaps1+=gaplen		
	
	#watch out: i have a feeling that if this fails occasionally,
	#it would be because I only check for sequence1's insertions and deletions
	#but I don't check for sequence2's insertions and deletions

	#TODO: construct a test where sequence2 also has some indels and see how this behaves.

	#lastly, we format the aligned sequences as FASTA (60 characters per line)
	fasta1.write(">hg19|chr21|produced using variants from "+member[2]+"|first sequence\n")
	fasta2.write(">hg19|chr21|produced using variants from "+member[2]+"|second sequence\n")
	fasta1.write("\n".join("".join(gappedsequence1[i:i+60]) for i in xrange(0, len(sequence1), 60))+"\n")
	fasta2.write("\n".join("".join(gappedsequence2[i:i+60]) for i in xrange(0, len(sequence2), 60))+"\n")
	
	vcffile.close()
	ref.close()
	window.close()
	
	fasta1.close()
	fasta2.close()
	
try:	
	configuration = sys.argv[1]
	mother = prepMember("Mother", configuration)
	father = prepMember("Father", configuration)
	child = prepMember("Child", configuration)
	VCFtoFASTA(child)
	VCFtoFASTA(mother)
	VCFtoFASTA(father)
except IndexError:
	print "Please select a configuration file."
except:
	print "The configuration file you selected is invalid."
	




