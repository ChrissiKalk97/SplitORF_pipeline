import sys

file = open(sys.argv[1],'r')
validlist={}
for line in file:
    elems = line.strip().split("\t")
    IDs=elems[4].split(",")
    for i in IDs:
        validlist[i]="Valid"
   

file2 = open(sys.argv[2],'r')
with open(sys.argv[3],'w') as f:
    for line in file2:
        elems = line.split("\t")
        foo = ([pos for pos, char in enumerate(elems[0]) if char == ":"])
        ID = elems[0][foo[0]+1:foo[1]]
        if ID in validlist.keys():
            f.write(line)
