#!/software/bin/python

# Author: lh3, converted to python and modified to add -C option by Aylwyn Scally
#
# About:
#   varfilter.py is a port of Heng's samtools.pl varFilter script into 
#   python, with an additional -C INT option. This option sets a minimum 
#   consensus score, above which the script will output a pileup line 
#   wherever it _could have_ called a variant, even if none is actually 
#   called (i.e. hom-ref positions). This is important if you want to
#   subsequently merge the calls with those for another individual to get a
#   synoptic view of calls at each site. Without this option, and in all 
#   other respects, it behaves like samtools.pl varFilter.
#   
#   Aylwyn Scally as6@sanger.ac.uk


# Filtration code:
#
# C low CNS quality (hom-ref only)
# d low depth
# D high depth
# W too many SNPs in a window (SNP only)
# G close to a high-quality indel (SNP only)
# Q low RMS mapping quality (SNP only)
# g close to another indel with higher quality (indel only)
# s low SNP quality (SNP only)
# i low indel quality (indel only)


import sys
import getopt

def usage():
	print '''usage: varfilter.py [options] [cns-pileup]

Options: -Q INT	minimum RMS mapping quality for SNPs
		 -q INT	minimum RMS mapping quality for gaps
		 -d INT	minimum read depth 
		 -D INT	maximum read depth
		 -S INT	minimum SNP quality
		 -i INT	minimum indel quality
		 -C INT	minimum consensus quality for hom-ref sites

		 -G INT	min indel score for nearby SNP filtering
		 -w INT	SNP within INT bp around a gap to be filtered

		 -W INT	window size for filtering dense SNPs
		 -N INT	max number of SNPs in a window

		 -l INT	window size for filtering adjacent gaps

		 -p print filtered variants'''

def varFilter_aux(first, is_print):
	try:
		if first[1] == 0:
			sys.stdout.write("\t".join(first[4:]) + "\n")
		elif is_print:
			sys.stderr.write("\t".join(["UQdDWGgsiCX"[first[1]]] + first[4:]) + "\n")
	except IOError:
		sys.exit()
 
mindepth = 3
maxdepth = 100
gapgapwin = 30
minsnpmapq = 25
mingapmapq = 10
minindelscore = 25
scorefactor = 100
snpgapwin = 10
densesnpwin = 10
densesnps = 2
printfilt = False
minsnpq = 0
minindelq = 0
mincnsq = 0

try:
	options, args = getopt.gnu_getopt(sys.argv[1:], 'pq:d:D:l:Q:w:W:N:G:S:i:C:', [])
except getopt.GetoptError:
	usage()
	sys.exit(2)
for (oflag, oarg) in options:
	if oflag == '-C':
		mincnsq = int(oarg)
	elif oflag == '-D':
		maxdepth = int(oarg)
	elif oflag == '-G':
		minindelscore = int(oarg)
	elif oflag == '-N':
		densesnps = int(oarg)
	elif oflag == '-Q':
		minsnpmapq = int(oarg)
	elif oflag == '-S':
		minsnpq = int(oarg)
	elif oflag == '-W':
		densesnpwin = int(oarg)
	elif oflag == '-d':
		mindepth = int(oarg)
	elif oflag == '-i':
		minindelq = int(oarg)
	elif oflag == '-l':
		gapgapwin = int(oarg)
	elif oflag == '-p':
		printfilt = True
	elif oflag == '-q':
		mingapmapq = int(oarg)
	elif oflag == '-s':
		scorefactor = int(oarg)
	elif oflag == '-w':
		snpgapwin = int(oarg)
inp = sys.stdin if len(args) < 1 else open(args[0])
# calculate the window size
max_dist = max(gapgapwin, snpgapwin, densesnpwin)

staging = []
for t in (line.strip().split() for line in inp):
	(flt, score) = (0, -1)
	# non-var sites
	if t[3] == '*/*':
		continue
	is_snp = t[2].upper() != t[3].upper()
	if not is_snp and not mincnsq:
		continue
	# clear the out-of-range elements
	while staging:
		# Still on the same chromosome and the first element's window still affects this position?  
		if staging[0][4] == t[0] and int(staging[0][5]) + staging[0][2] + max_dist >= int(t[1]):
			break
		varFilter_aux(staging.pop(0), printfilt)

	# first a simple filter
	if int(t[7]) < mindepth:
		flt = 2
	elif int(t[7]) > maxdepth:
		flt = 3
	if t[2] == '*': # an indel
		if minindelq and minindelq > int(t[5]):
			flt = 8
	elif is_snp:
		if minsnpq and minsnpq> int(t[5]):
			flt = 7
	elif mincnsq > int(t[4]):
		flt = 9

	# site dependent filters
	dlen = 0
	if flt == 0:
		if t[2] == '*':
			# If deletion, remember the length of the deletion
			(a,b) = t[3].split('/')
			alen = len(a) - 1
			blen = len(b) - 1
			if alen>blen:
				if a[0] == '-': dlen=alen 
			elif b[0] == '-': dlen=blen 

			if int(t[6]) < mingapmapq:
				flt = 1
			# filtering SNPs
			if int(t[5]) >= minindelscore:
				for x in (y for y in staging if y[3]):
					# Is it a SNP and is it outside the SNP filter window?
					if x[0] >= 0 or int(x[5]) + x[2] + snpgapwin < int(t[1]):
						continue
					if x[1] == 0:
						x[1] = 5

			# calculate the filtering score (different from indel quality)
			score = int(t[5])
			if t[8] != '*':
				score += scorefactor * int(t[10])
			if t[9] != '*':
				score += scorefactor * int(t[11])
			# check the staging list for indel filtering
			for x in (y for y in staging if y[3]):
			  # Is it a SNP and is it outside the gap filter window
				if x[0] < 0 or int(x[5]) + x[2] + gapgapwin < int(t[1]):
					continue
				if x[0] < score:
					x[1] = 6
				else:
					flt = 6
					break
		else:
			if int(t[6]) < minsnpmapq:
				flt = 1
			k = 1 + sum(
				1
				for x in (y for y in staging if y[3])
				if x[0] < 0
				and int(x[5]) + x[2] + densesnpwin >= int(t[1])
				and x[1] in [0, 4, 5]
			)
			# filtering is necessary
			if k > densesnps:
				flt = 4
				for x in (y for y in staging if y[3]):
					if x[0] < 0 and int(x[5]) + x[2] + densesnpwin >= int(t[1]) and x[1] == 0:
						x[1] = 4
			else: # then check gap filter
				for x in (y for y in staging if y[3]):
					if x[0] < 0 or int(x[5]) + x[2] + snpgapwin < int(t[1]):
						continue
					if x[0] >= minindelscore:
						flt = 5
						break

	staging.append([score, flt, dlen, is_snp] + t)

# output the last few elements in the staging list
while staging:
	varFilter_aux(staging.pop(0), printfilt)
