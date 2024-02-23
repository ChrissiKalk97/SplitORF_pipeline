# This file is used to merge the bedfile entries of the unique region bedfiles in order to
# correct for the bias introduced through the alignment with MUMmer. Two succeding entries
# will be merged together if the end position of the first and the starting position of the second
# are less than X apart. Where X is a parameter set by the user and should be set at the same
# value as the minimum length parameter of the MUMmer maxmatch call.

import sys

# Bedfile to be merged
infile = open(sys.argv[1], "r")
lines = infile.readlines()
# Parameter X ()should be set at the same value as the minimum length parameter of the MUMmer maxmatch call.
gap = int(sys.argv[2])
previous = 3*[""]
# Resulting bedfile
with open(sys.argv[3], "w") as out:
    for i in lines:
        columns = i.split()
        if columns[0]==previous[0]:
            if int(columns[1])-int(previous[2]) <= gap:
                previous[2]=columns[2]
            else:
                if (int(previous[2])-int(previous[1]) >= gap):
                    out.write(previous[0] + "\t" + previous[1] + "\t" + previous[2] + "\n")
                previous = columns
        else:
            if (previous[0]!=""):
                if (int(previous[2])-int(previous[1]) >= gap):
                    out.write(previous[0] + "\t" + previous[1] + "\t" + previous[2] + "\n")
            previous=columns
    if (int(previous[2])-int(previous[1]) >= gap):
        out.write(previous[0] + "\t" + previous[1] + "\t" + previous[2] + "\n")