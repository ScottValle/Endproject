gile = open("contig_clustering_all.csv", "r")
file = open("genomes.tsv", "r")
tile = open("genomesspecies.tsv")
oile = open("genomesnamed.tsv", "w")
boile = open("genomesspeciesnamed.tsv", "w")
#opening all the files
anidict = {}
keydict = {}
for line in gile.readlines()[1:]:
    line =line.split(',')
    if anidict.get(line[23]):
        anidict[line[23]]["host"].add(line[5].replace("_", " ")), anidict[line[23]]["genus"].add(line[17].rstrip('"') + ' ' + line[18].lstrip('"'))
    else:
        anidict[line[23]] = {"host": {line[5].replace("_", " ")}, "genus": {line[17].rstrip('"') + ' ' + line[18].lstrip('"')}}
    keydict[line[0].strip('"')] = line[23]
#collects all the combinations of host genus combinations from the species
text = file.readlines()
oile.write("\t".join(text[0].strip("\n").split("\t")) + '\t"Genus"\t"Host"' + '\n')
for line in text[1:]:
    line = line.split("\t")
    line[-1] = line[-1].strip("\n")
    if list(anidict[keydict.get(line[0].strip('"'))]["genus"])[0] == "NA NA":
        line.append(list(anidict[keydict.get(line[0].strip('"'))]["genus"])[1].replace('"', ""))
    else:
        line.append(list(anidict[keydict.get(line[0].strip('"'))]["genus"])[0].replace('"', ""))
    temp = list(anidict[keydict.get(line[0].strip('"'))]["host"])
    if "NA" in temp:
        temp.remove("NA")
        if len(temp) == 0:
            temp.append('"Unknown"')
    if '"Unknown"' in temp and len(temp) > 1:
        temp.remove('"Unknown"')
    line.append("_".join(temp).replace('"', ""))
    #print("\t".join(line))
    oile.write("\t".join(line) + '\n')
#writes the annotated and named files which could have been done with a loop, but this really needs to be run just the once, so this works
rext = tile.readlines()
boile.write("\t".join(rext[0].strip("\n").split("\t")) + '\t"Genus"\t"Host"' + '\n')
for line in rext[1:]:
    line = line.split("\t")
    line[-1] = line[-1].strip("\n")
    if list(anidict[keydict.get(line[0].strip('"'))]["genus"])[0] == "NA NA":
        line.append(list(anidict[keydict.get(line[0].strip('"'))]["genus"])[1].replace('"', ""))
    else:
        line.append(list(anidict[keydict.get(line[0].strip('"'))]["genus"])[0].replace('"', ""))
    temp = list(anidict[keydict.get(line[0].strip('"'))]["host"])
    if "NA" in temp:
        temp.remove("NA")
        if len(temp) == 0:
            temp.append('"Unknown"')
    if '"Unknown"' in temp and len(temp) > 1:
        temp.remove('"Unknown"')
    line.append("_".join(temp).replace('"', ""))
    #print("\t".join(line))
    boile.write("\t".join(line) + '\n')



