
# given a flycatcher scaffold map, this will convert swallow scaffolds and position in a vcf to flycatcher positions
# this involves substantial sorting of potentially very large amounts of data. Good to run in separate directory, because of temp files.

# NOTE: do not redirect the output, because this file automatically writes to an output file

# NOTE: need to run "ulimit -n 10000" first on some sysmtems

# input
# 1. flycatcher map
# 2. the vcf

import sys, os, math
from operator import itemgetter


# first read the map
mapp = {} # dictionary of swallow scaffolds as keys with associated flycatcher hit info as values
lgMap = {} # numbering flycatcher scaffolds, because the previous scaffold IDs are too long
scafCounter = 1
scaffoldLengths = {} # initializing a dictionary of scaffold lengths
with open(sys.argv[1], "r") as infile:
    for line in infile:
        newline = line.strip().split()
        swalScaf, flyScaf, flyPos, chrom = newline[:]
        mapp[swalScaf] = [flyScaf, flyPos, chrom]
        if flyScaf not in lgMap: # counting the number of flycatcher scaffolds on each chromosome
            lgMap[flyScaf] = scafCounter
            scaffoldLengths[str(scafCounter)] = 0 # collecting the lengths later
            scafCounter += 1
            
sys.stderr.write("\nmap read, reading vcf\n")

# read in vcf, and replace swallow positions with flycatcher positions
flyPosVcf = {}
tempFileName = "unsorted.tmp"
tempFile = open(tempFileName, "w")
with open(sys.argv[2], "r") as infile:
    for line in infile:
        if line[0] == "#":
            header = line
        else:
            newline = line.strip().split()
            scaf, pos = newline[0:2]
            if scaf in mapp: # if the swallow scaf has a flycatcher hit
                scafID = mapp[scaf][0]
                scafNum = str(lgMap[scafID])
                #
                outline = list(newline)
                outline[0] = scafNum # replace the chrom column with flycatcher scaffold number
                newPos = int(mapp[scaf][1]) + int(pos) - 1 # replace the pos with the flycatcher pos, + the swallow bp
                outline[1] = newPos 
                tempFile.write('\t'.join(map(str,outline)) + "\n") # write to unsorted temp file
                #
                if newPos > scaffoldLengths[scafNum]:
                    scaffoldLengths[scafNum] = newPos        
tempFile.close() # close the temp file after writing is complete

sys.stderr.write("\nvcf read. unsorted temp file written, now re-reading the mapped temp file and splitting into smaller temp files \n")


# now open up a bunch of temporary files, apprx 1 per megaBase on each scaffold = 1000 files, and write data to them
tempFiles = {} # big dictionary of tempFile paths
for scaf in lgMap:
    scafNum = str(lgMap[scaf])
    scafLength = scaffoldLengths[scafNum]
    numWindows = int(math.ceil(float(scafLength) / 1000000))
    tempFiles[scafNum] = {}
    for w in range(numWindows):
        windowNum = w+1
        tempFiles[scafNum][str(windowNum)] = open("temp_scaf" + scafNum + "_" + str(windowNum) + ".tmp", "w")    
with open(tempFileName, ) as infile: # now read in the mapped, unsorted vcf, writing to individual temporary files for each scaffold
    for line in infile:
        newline = line.strip().split()
        scaf = newline[0]
        pos = int(newline[1])
        windowNum = (pos / 1000000) + 1 # dividing integers gives interger, rounded down
        tempFiles[scaf][str(windowNum)].write(line)
for scaf in tempFiles: # close all temp files, done with writing to them
    for tfile in tempFiles[scaf]:
        tempFiles[scaf][tfile].close()
bashCommand = "rm "+tempFileName # remove the big temp file
os.system(bashCommand)

sys.stderr.write("\ntemp files written. now re-reading the individual temp files and sorting them, adding to final VCF\n")
        
# read in individual temp files, sort, rewrite to final vcf
origVcfName = sys.argv[2]
outName = ".".join(origVcfName.split(".")[:-1]) + "_flycatcherScaf.vcf"
outputVcf = open(outName, "w")
outputVcf.write(header)
for scafNum in range(1, len(tempFiles)+1):
    scaf = str(scafNum)
    for windowNum in range(1, len(tempFiles[scaf])+1):
        window = str(windowNum)
        tempFname = "temp_scaf" + scaf + "_" + window + ".tmp"
        unsortedData = []
        with open(tempFname, "r") as infile:
            for line in infile:
                newline = line.strip().split()
                newline[1] = int(newline[1])
                unsortedData.append(newline)
        orderedList = sorted(unsortedData, key=itemgetter(1)) # sort
        for snp in orderedList:
            outputVcf.write("\t".join(map(str, snp)) + "\n")
        bashCommand = "rm "+tempFname # remove the smaller temp files
        os.system(bashCommand) 
outputVcf.close()







