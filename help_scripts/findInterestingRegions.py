from collections import deque
#import time

child = open("NA12880-clean.vcf", "r")
window = open("NA12880-clean-INDELwindows.vcf", "w")

win = deque([])

maxwin = 0
where = 0,'0'
#winnum = [0,0,0,0,0,0,0]
prevchr = '0'
for line in child:
	if line[0]=='#':
		window.write(line)
		continue
	if line[0]=='X':
		break
	info = line.split()
	pos = int(info[1])
	if len(info[3])==len(info[4]):
		continue
	if info[9][0:3]=="1/1": #homozygous, doesn't count
		continue
	if not info[0]==prevchr:
		win = deque([])
	#print len(win),
	tags = info[8]
	vals = info[9]
	win.append(pos)
	uppest = win.popleft()
	while pos - uppest > 1000 and len(win):
		uppest = win.popleft()
	win.appendleft(uppest)
	windowsize = len(win)
	if windowsize>=8 and info[0]=='21':
		print windowsize, (pos,info[0])
	if windowsize > maxwin:
		maxwin = windowsize
		where = (pos, info[0])
		#print maxwin, where
	#winnum[windowsize-1] = winnum[windowsize-1]+1
	#print str(windowsize)
	#print vals+":"+str(windowsize)
	#time.sleep(2)
	tags = tags + ":WI"
	vals = vals + ":"+str(windowsize)
	info[8] = tags
	info[9] = vals
	line = '\t'.join(info)+"\n"
	window.write(line)
	prevchr = info[0]

print "The biggest number of variants in any 1000 consecutive positions is " + str(maxwin) +" at "+str(where)

#for index in xrange(7):
#	print "There are "+str(winnum[index])+" windows of variants that contain exactly "+ str(index+1)+ " variants"

child.close()
window.close()
