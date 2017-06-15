# python script.py > dayhoff.fig80.txt

fn = 'dayhoff.corrected.changes.txt'
FH = open(fn)

blocks = FH.read().strip().split('\n\n')
data = blocks[0].strip().split('\n')

pL = list()
aaL = list()

for line in data:
    L = line.strip().split()
    aa = L.pop(0)
    aaL.append(aa)
    print aa, ' ',
    L = [e.rjust(3) for e in L]
    print '  '.join(L)
    
print '      ' + '    '.join(aaL[:-1])