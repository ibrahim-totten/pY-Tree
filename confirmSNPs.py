import csv

def genotypeAtCluster(rsid):   
    f = open('gen.csv', 'r')
    reader = csv.reader(f)
    #Search and return genotype             Add binary search option for sorted rsid csvs
    for row in reader:
        if row[0] == rsid:
            f.close()
            return row[3]
    f.close()
    return None

def confirmGenotype(genotype,rsid):
    return genotypeAtCluster(rsid) == genotype

def SNPtoRSID(snp):
    f = open('SNPIndex.csv', 'r')
    reader = csv.reader(f)
    #Search and return rsid             Add binary search option for sorted rsid csvs
    for row in reader:
        if row[0] == snp:
            f.close()
            return row[3]
    f.close()
    return None

def confirmSNP(snp):
    f = open('SNPIndex.csv', 'r')
    reader = csv.reader(f)
    #Search and return rsid             Add binary search option for sorted rsid csvs
    for row in reader:
        if row[0] == snp:
            f.close()
            return confirmGenotype(getDeriv(row[6]), row[3])
    f.close()
    return None

def confirmRSID(rsid):
    pass

def getDeriv(mutation):
    if len(mutation) == 4:
        return mutation[3]

def getAliases(snp):
    f = open('SNPIndex.csv', 'r')
    reader = csv.reader(f)
    #Search and return rsid             Add binary search option for sorted rsid csvs
    for row in reader:
        if row[0] == snp:
            f.close()
            return row[2].split("; ")
    f.close()
    return None

def getRsidsFromBranch(branchName):
    f = open('SNPIndex.csv', 'r')
    reader = csv.reader(f)
    #Search and return rsid             Add binary search option for sorted rsid csvs
    rsids = []
    alreadyChecked = []
    
    for row in reader:
        if row[1] == branchName and row[3] not in alreadyChecked:
            rsids.append((row[3],row[6]))
            alreadyChecked.append(row[3])
    f.close()
    return rsids


#print(genotypeAtCluster("rs774860428"))
#print(confirmGenotype( "T", "rs774860428"))
print(getRsidsFromBranch("R"))
print(getAliases("BT"))
print(confirmSNP("L146"))