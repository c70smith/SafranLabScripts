

# This is doing a calculation according to the PLtoGP function in the VariantAnnotation package (BioConductor)
# https://www.rdocumentation.org/packages/VariantAnnotation/versions/1.18.5/topics/GLtoGP
# GLtoGP computes the probability of each genotype as 10^x / sum(10^x). PLtoGP first divides by -10 and then proceeds as in GLtoGP.
# i.e. this calculation is 10^(PL/-10) / sumForEachOfTheGenotypes( 10^(PL/-10) ), then
#      converting to a single value between 0 and 2 (what's called the "mean" in Alex's script) by adding 1*the hetero GP plus 2* the ALT GP
# The reason I'm writing my own code is because the calculation is quite simple and this is way easier to use, because the R package requires certain metadata-header lines explaining the fields, and it will take coding to convert the three* GPs output by the R package to one value anyways.
# And this is the exact same as Alex's script, but bypasses the intermediate files and programs, and we can skip referencing him.


import sys, numpy

with open(sys.argv[1], 'r') as vcf:

    # print header
    header = vcf.readline().strip().split()
    newHead = header[0:2] + header[9:]
    print "\t".join(newHead)
    
    # go through vcf and convert PLs to GPs
    for line in vcf:
        newline = line.strip().split()
        samples = newline[0:9]
        genos = newline[9:]
        genoScores = [None]*len(genos)
        for geno in range(len(genos)):
            data = genos[geno].split(":")
            converted = [None,None,None]
            pls = data[1].split(",")
            if len(pls) == 3:
                for pl in range(3):
                    converted[pl] = 10**(float(pls[pl]) / -10)
                if converted[0] == converted[1] and converted[0] == converted[2]: # all equally likely, output NA
                    GP = "NA"
                else:
                    GP =  round( (converted[1] / sum(converted)) + 2*(converted[2]/sum(converted)) , 5)
                genoScores[geno] = GP
            else:
                sys.stderr.write("\n\nbiallelic only\n\n")
                1/0
        outline = samples[0:2] + genoScores
        print "\t".join(map(str, outline))
            
                    
            




