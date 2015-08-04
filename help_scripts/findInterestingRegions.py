from collections import deque
import ConfigParser
import sys

def prepVCF(configpath):
	config = ConfigParser.RawConfigParser()
	config.read(configpath)
	child = config.get("Child", "vcf_path")
	thres = int(config.get("Threshold", "value"))
	res = config.get("Regions", "path")
	return [child, thres, res]

def getInterestingRegions (childPath, threshold, regionPath):
	child = open(childPath, "r")
	region = open(regionPath, "w")
	win = deque([])
	prevchr = '0'
	for line in child:
		if line[0]=='#': #a header, no variant information
			continue
		info = line.split()
		pos = int(info[1])
		if len(info[3])==len(info[4]): #a replacement, not an INDEL, doesn't count
			continue
		if info[9][0:3]=="1/1": #homozygous, doesn't count
			continue
		if not info[0]==prevchr:
			win = deque([])
		tags = info[8]
		vals = info[9]
		win.append(pos)
		uppest = win.popleft()
		while pos - uppest > 1000 and len(win):
			uppest = win.popleft()
		win.appendleft(uppest)
		windowsize = len(win)
		if windowsize>=threshold:
			region.write( "----- Region with " + str(windowsize) + " interesting variants: -----\n")
			region.write(info[0]+"\n")
			region.write(str(pos-1000)+"\n")
			region.write(str(pos)+"\n")
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
#except:
#	print "The configuration file you selected is invalid."
	
