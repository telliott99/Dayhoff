import sys
FH = open('data/dayhoff.fig80.txt', 'r')
data = FH.read().strip().split('\n')
FH.close()

#--------------------------------------------

# for convenience, save aa in order as aaL
# helps when grabbing the correct column, below

aaL = data.pop().strip().split()

# keys for sorting
def f0(s):   return aaL.index(s[0])
def f1(s):   return aaL.index(s[1])

# write aa order to a file for help with R code
def save_order():
    FH = open('dayhoff.aaorder.txt', 'w')
    FH.write(' '.join(aaL) + '\n' )
    FH.close()

# save_order()

#--------------------------------------------

# Dayhoff data is a triangular table of totalChanges
# save in a dict
# note changeD[k1k2] = changeD[k2k1]

changeD = dict()
for line in data:
    L = line.split()
    aa1 = L.pop(0)
    # print aa1, L
    for i, value in enumerate(L):
        aa2 = aaL[i]
        n = int(value)
        changeD[aa1+aa2] = n
        changeD[aa2+aa1] = n

# total totalChanges per aa in same order as aaL
totalChanges = list()
for aa in aaL:
    L = [k for k in changeD if aa == k[0]]
    totalChanges.append(
        sum([changeD[k] for k in L]))
    
def show():
    for i,aa in enumerate(aaL):
        # all keys w/aa as first letter of key
        L = [k for k in changeD if aa == k[0]]
        # sort by pos of 2nd aa in aaL
        L = sorted(L, key=f1)
        # get the actual values
        L = [changeD[k] for k in L]
        print aa,totalChanges[i],sum(L), L

# matches the paper and entered data
# show();  sys.exit()

#--------------------------------------------

# aa frequency data from the paper
# in order from most to least frequent

FH = open('data/dayhoff.frequencies.txt', 'r')
data = FH.read().split('\n')
FH.close()

# read from data in same order as aaL
def getFreq(aa):
    for line in data:
        if aa == line[0]:
            return float(line.split()[1])
freq = [getFreq(aa) for aa in aaL]

def write_freq():
    # write to a file for help with R code
    FH = open('dayhoff.frequencies.v2.txt', 'w')
    L = [str(freq[i]) for i in range(20)]
    FH.write('\t'.join(L) + '\n')
    FH.close()

# write_freq()
 
#--------------------------------------------

# mutability is totalChanges observed / frequency
# normalized to alanine = 100

L = [totalChanges[i]*1.0/freq[i] for i in range(20)]
mutAlanine = L[0]
mutability = [e * 100 / mutAlanine for e in L]

def show2():   
    print sum(totalChanges), 'totalChanges'
    print sum(changeD.values()), 'changeD values'
    for i in range(20):
        print aaL[i], str(totalChanges[i]).rjust(6), 
        print '%3.3f' % freq[i],
        print str(int(round(mutability[i]))).rjust(4)

# show2();  sys.exit()
 
#--------------------------------------------

# we've shown we can get close to their values
# for "mutability", now we'll actually use theirs

FH = open('data/dayhoff.mutabilities.txt', 'r')
data = FH.read().split('\n')
FH.close()
L = [e.split('\t') for e in data]
L = [(e[0],int(e[1])) for e in L]
D = dict(L)
L = [D[aa] for aa in aaL]
# print L
 
#===========================================+

# calculate PAM1
# opposite of current convention:
# cols indexed by j for current
# rows indexed by i for future

# go back to the change dictionary cD
PAM1 = dict()

for oldAA in aaL:
    total = 0
    for k in changeD:
         # order doesn't matter here, pick k[0]
        if not k[0] == oldAA:   continue
        newAA = k[1]
        nchanges = changeD[k]
        
        i = aaL.index(newAA)
        j = aaL.index(oldAA)
        S = totalChanges[j]  # all changes involving oldAA
        mut = mutability[j]  # mutability of oldAA
        
        K = 1.33
        value = K * mut * nchanges / S
        # print k, nchanges, S, mut,'%3.2f' % value
        # print '%3.2f' % ((K * nchanges) / (418.85 * freq[j]))
        total += value
        PAM1[k] = int(round(value))
    PAM1[oldAA+oldAA] = 10000 - int(total)


print ' ' + ''.join([aa.rjust(5) for aa in aaL])
print
for newaa in aaL:
    print newaa,
    for oldaa in aaL:
        print str(PAM1[oldaa+newaa]).rjust(4),
    print
    print
         
#===========================================+
    
# print PAM to a file in their order of aa:
FH = open('PAM1.txt','w')
FH.write('\t'.join(aaL) + '\n')
lines = list()
for aa2 in aaL:
	line = list()
	for aa1 in aaL:
		line.append( ('%i' % PAM1[aa1+aa2]).rjust(6) )
	lines.append(''.join(line))
FH.write('\n'.join(lines))
   
    



