import sys
file = open(sys.argv[1],'r')
with open(sys.argv[2], 'w') as f:
    for line in file:
        elems = line.strip()
        elems = elems.split("\t")
        if int(elems[6]) == 1:
            start = int(elems[2]) - int(elems[4])
            end = int(elems[3]) - int(elems[4])
            f.write(elems[0] + "|"+elems[1] + "\t" + str(start) + "\t" + str(end) + "\n")
        else:
            start = int(elems[5]) - int(elems[3])
            end = int(elems[5]) - int(elems[2])
            f.write(elems[0] + "|" + elems[1] + "\t" + str(start) + "\t" + str(end) + "\n")