from collections import deque
import ConfigParser
import sys

def prepVCF(configpath):
	"""
	Assuming configpath is a config file with sections Child, Threshold and Regions.
	Child's vcf_path describes the VCF file of the child of the trio. Regions's path
	describes the path to the output file. Threshold's value is an integer that
	describes what is the minimum number of "interesting" variants that has to be
	present in the 1000bps window so that this window makes it to the output Results file.	
	"""
	config = ConfigParser.RawConfigParser()
	config.read(configpath)
	child = config.get("Child", "vcf_path")
	thres = int(config.get("Threshold", "value"))
	res = config.get("Regions", "path")
	return [child, thres, res]

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
		not line.split()[4] == "."
	):
		return True
	else:
		return False
	

def getInterestingRegions (childPath, threshold, regionPath):
	"""
	Assumming childPath is a path to a VCF file, threshold is an integer and regionPath
	is a path to a text file. We look through the VCF file and find the "interesting"
	variants. An "interesting" variant is a variant that is a heterozygous INDEL
	(insertion or deletion). These "interesting" variants are difficult to phase,
	since homozygous variants, naturally, occur in both the mother and the father
	copy of the chromosome, and SNP replacements of one character (e.g. A->C) can be
	phased with tools like WhatsHap.
	
	The resulting file, recorded in regionPath, gives information about the number of the
	interesting variants in a window of 1000 bps and then presents the chromosome,
	the start and the end of the window in the format of the .posinfo file. The idea is
	that one can copy a window from the regionPath file straight to the .posinfo file.
	"""
	child = open(childPath, "r")
	region = open(regionPath, "w")
	win = deque([]) #save the interesting variant positions in radius of 1000 bps
	prevchr = '0' #the chromosome of the previous line
	for line in child:
		if line[0]=='#': #a header, no variant information
			continue
		if not doesLineCount(line):
			continue
		info = line.split()
		pos = int(info[1])
		if len(info[3])==len(info[4]): #a replacement, not an INDEL, doesn't count
			continue
		if info[9][0:3]=="1/1": #homozygous, doesn't count
			continue
		if not info[0]==prevchr: #if you get on a new chromosome, start counting over
			win = deque([])
		win.append(pos)
		uppest = win.popleft()
		while pos - uppest > 1000 and len(win): #update how many interesting variants are there in radius of 1000 bps
			uppest = win.popleft()
		win.appendleft(uppest)
		windowsize = len(win) #how many interesting (heterozygous INDELs) variants exist in these 1000 bps
		if windowsize>=threshold:
			region.write( "----- Region with " + str(windowsize) + " interesting variants: -----\n")
			region.write(info[0]+"\n")
			region.write(str(pos-1000)+"\n")
			region.write(str(pos+1)+"\n")
		prevchr = info[0]

	child.close()
	region.close()

try:	
	configuration = sys.argv[1]
	setting = prepVCF(configuration)
	try:
		getInterestingRegions(setting[0], setting[1], setting[2])
	except IOError:
		print "Please make sure the path to the resulting FASTA file is valid."
except IndexError:
	print "Please select a configuration file."
except:
	print "The configuration file you selected is invalid."
	raise
	
