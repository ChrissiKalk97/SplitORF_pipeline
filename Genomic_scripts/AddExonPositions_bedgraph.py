import sys
file = open(sys.argv[1],'r')
dict={}
for line in file:
    elems = line.strip()
    elems = elems.split("\t")
    ids=elems[0]
    temp=[int(elems[1]),int(elems[2])]
    if ids in dict.keys():
        t=dict[ids]
        t.append(temp)
        dict[ids]=t
    else:
        dict[ids]=[temp]
for k, v in dict.items():
        entries = sorted(v)
        end=entries[0][1]-entries[0][0]+1
        add=entries[0][0]
        entries[0]=[end,add]
        last=end
        for i in range(1, len(entries)):
            interval = entries[i][1]-entries[i][0]+1
            end=last+interval+1
            add = entries[i][0]-(last+1)
            entries[i]=[end,add]
            last=end
        dict[k]= entries

file = open(sys.argv[2],'r')
dict2={}
for line in file:
    elements = line.strip()
    elements = elements.split("\t")
    ids2=elements[0]
    #print(elements[0] + "\t" + elements[1] + "\t" + elements[2] + "\t" + elements[3])
    temp2=[int(elements[1]),int(elements[2]),int(elements[3])]
    if ids2 in dict2.keys():
        t2=dict2[ids2]
        t2.append(temp2)
        dict2[ids2]=t2
    else:
        dict2[ids2]=[temp2]
for k, v in dict2.items():
    #print(k)
#    if k == "ENSG00000124193|ENST00000483871":
    for i in v:
        j = 0
        while i[0] > dict[k][j][0]:
            j=j+1
        i[0]=i[0]+dict[k][j][1]
        j = 0
        while i[1] > dict[k][j][0]:
            #print(i[1], dict[k])
            j=j + 1
        i[1] = i[1] + dict[k][j][1]
with open(sys.argv[3], 'w') as f:
    for k, v in dict2.items():
        for i in v:
            f.write(k + "\t" + str(i[0]) + "\t" + str(i[1]) + "\t" + str(i[2]) + "\n")
#test=[106,514]
#test2=[761,1409]
#list=[]
#list.append(test)
#list.append(test2)
#addlist=[[212,0],[361,220],[630,863],[755,1212],[964,1349],[1048,1622],[1680,1726]]
#for i in list:
#    j = 0
#    while i[0] > addlist[j][0]:
#        j=j+1
#    i[0]=i[0]+addlist[j][1]
#    j = 0
#    while i[1] > addlist[j][0]:
#        j=j + 1
#    i[1] = i[1] + addlist[j][1]
#print(list)