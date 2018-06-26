

import sys
with open(sys.argv[1], "r") as infile:
    header = infile.readline()
    for line in infile:
        newline = line.strip().split()
        col1 = "-".join(newline[0:2])
        print "\t".join( [col1] + newline[2:] )

        
