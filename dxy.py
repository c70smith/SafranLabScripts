


# NOTE:
# the input vcf should be unfiltered, but should contain only the samples you want to analyze
# currently calculating dxy for windows with as few as one snp
# the dxy calculation requires haplotypes, which we do not have, therefore I'm counting the total allelic differences between two haploids 
# also we are ignoring snps without sufficient coverage which could make things weird, or downwardly bias dxy

# input
# 1. vcf
# 2. population / subspecies map, with sample ID (corresponding to those in the VCF), and popuation ID
# 3. min reads to genotype a sample
# 4. min proportion of samples covered from each population to calcualte dxy
# 5. the number of columns in the vcf before the genotype columns begin (usually 9)
# 6. window size
# 7. step size

import sys, numpy
minReadDepth = int(sys.argv[3])
minPropSamps = float(sys.argv[4])
numInfoCols = int(sys.argv[5])
wSize = int(sys.argv[6])  
sSize = int(sys.argv[7])  






# function for calculating Dxy from a set of snps
def processWindow(theData, populDict, sampleOrder, reverseDict):

    # make a dictionary of haplotyes, with their frequencies, for each population
    haplos1 = {} 
    haplos2 = {}
    for sampl in range(len(theData)):
        popul = populDict[sampleOrder[sampl]]
        hap = ''.join(map(str,theData[sampl]))
        if popul == 0:
            if hap not in haplos1:
                haplos1[hap] = 0
            haplos1[hap] += 1
        elif popul == 1:
            if hap not in haplos2:
                haplos2[hap] = 0
            haplos2[hap] += 1
        else:
            sys.stderr.write("\n\nHouston we have a problem\n\n")
            1/0

    # compare each pair of haplotypes
    dxy = 0
    for s1 in haplos1:
        h1 = list(s1)
        if set(h1) != set(["N"]): # missing data only
            x = float(haplos1[s1]) / reverseDict[0]
            #a[:] = [4 if x==1 else x for x in a]
            #h1[:] = [x for x in h1 if x != "N"]
            for s2 in haplos2:
                h2 = list(s2)
                if set(h2) != set(["N"]):
                    y = float(haplos2[s2]) / reverseDict[1]
                    c1 = []
                    c2 = []
                    for p in range(len(h1)): # going through to remove Ns
                        if h1[p] == "N" or h2[p] == "N":
                            pass
                        else:
                            c1.append(h1[p])
                            c2.append(h2[p])
                            
                    diffs = sum(abs(numpy.subtract( map(int,c1), map(int,c2) )))
                    dxy += (diffs * x * y)
    
    return dxy
    
    







# function for processing the data from one whole chromosome
def processChrom(snps, windowSize, stepSize, popDict, minDepth, minProp, numInfo, sampleList, reversePopDict, chromo):
    finalPosition = int(snps[-1][1])
    numWindows = ((finalPosition - windowSize) / stepSize) + 1 # this should do it
    startPosition = -stepSize + 1 # snps are 1-indexed
    endPosition = windowSize - stepSize # initial value only
    snpsInsideWindow = []
    currentSnpIndex = 0
    numSnps = len(snps)
    goodSnps = []
    numSamples = len(sampleList)
    for window in range(numWindows):
        startPosition += stepSize 
        endPosition += stepSize

        # use while loops to prune snps no longer inside window, and add new snps
        if snpsInsideWindow != []: # if not empty, prune
            keepOnSliding = True
            while keepOnSliding == True:
                if int(snpsInsideWindow[0][1]) < startPosition:
                    snpsInsideWindow.pop(0)
                else:
                    keepOnSliding = False
                if snpsInsideWindow == []: 
                    keepOnSliding = False
        keepOnSliding = True
        while keepOnSliding == True:
            if currentSnpIndex <= (numSnps-1):
                if int(snps[currentSnpIndex][1]) <= endPosition:
                    snpsInsideWindow.append(snps[currentSnpIndex])
                    currentSnpIndex += 1
                else:
                    keepOnSliding = False
            else:
                keepOnSliding = False

        # apply additional checks, and construct haplotypes
        numSnpsInWindow = len(snpsInsideWindow)
        if numSnpsInWindow > 0:
            removeList = []
            snpMat = [ ["N"]*numSnpsInWindow for i in range(numSamples)] # m rows x n cols, m = num samples, n = num snps
            for snpPos in range(numSnpsInWindow):
                snp = snpsInsideWindow[snpPos]
                altAllele = snp[4]
                if altAllele in ["A","C","T","G"]:
                    data = snp[numInfo:]
                    popCounts = [0, 0]
                    for sample in range(numSamples):
                        splitSample = data[sample].split(":")
                        geno, dp = splitSample[0], int(splitSample[2])
                        if dp >= minDepth and geno != "./.":
                            p = popDict[header[sample]]
                            popCounts[p] += 1
                            alleles = sum(map(int, geno.split("/"))) # summing the vcf geno output to get 0, 1 or 2
                            snpMat[sample][snpPos] = alleles
                    for popu in range(2):
                        if (float(popCounts[popu]) / reversePopDict[popu]) < minProp:
                            removeList.append(snpPos)
            goodSnps = [ [] for i in range(numSamples)]
            for sample in range(numSamples):
                for snpPos in range(numSnpsInWindow):
                    if snpPos not in removeList:
                        goodSnps[sample].append(snpMat[sample][snpPos])
            if goodSnps[0] != []:
                d = processWindow(goodSnps, popDict, sampleList, reversePopDict)
                print "\t".join( map(str,[chromo, startPosition, endPosition, d]) )
    return









# main
with open(sys.argv[2], "r") as infile: # first, the pop map
    pDict = {}
    countPops = {}
    popNum = 0
    rDict = {}
    for line in infile:
        sampId, popId = line.strip().split()
        if popId not in countPops:
            countPops[popId] = popNum
            rDict[popNum] = 0
            popNum += 1
        pDict[sampId] = countPops[popId]
        rDict[countPops[popId]] += 1
    if len(countPops) != 2:
        sys.stderr.write("\n\n need exactly 2 populations \n\n")
        1/0

with open(sys.argv[1], "r") as infile: # now read vcf
    header = infile.readline().strip().split()
    if header[4] != "ALT":
        sys.stderr.write("\n\n fifth column needs to specify alternate allele \n\n")
        1/0
    header = header[numInfoCols:]
    
    # begin by reading in an entire chromosome at a time. I think this is the simplest way to code this, without reading in entire VCF which requires a lot of RAM
    newline = infile.readline().strip().split()
    currentChrom = newline[0]
    data = [newline]
    for line in infile:
        newline = line.strip().split()
        chrom = newline[0]
        if chrom == currentChrom:
            data.append(newline)
        else:
            results = processChrom(data, wSize, sSize, pDict, minReadDepth, minPropSamps, numInfoCols, header, rDict, currentChrom)
            currentChrom = chrom
            data = [newline]
    results = processChrom(data, wSize, sSize, pDict, minReadDepth, minPropSamps, numInfoCols, header, rDict, currentChrom) # processing the last chromosome outside the loop

