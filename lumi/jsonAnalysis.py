import os,itertools,sys

def incRange(a,b) : return range(a,b+1)
def lumis(pairs) : return set(sum([incRange(*pair) for pair in pairs],[]))
def pairs(lumis) : return [ (lambda l: [l[0][1],l[-1][1]] )(list(g))  for k,g in itertools.groupby( enumerate( sorted(lumis) ) , lambda (i,x):i-x ) ]

def intersect(a,b) : return dict( filter(lambda (k,p):p, [(k, pairs( lumis(a[k]) & lumis(b[k]) )) for k in set(a)&set(b) ] ) )
def union(a,b) : return dict( [ (k, pairs( lumis(a[k] if k in a else []) | lumis(b[k] if k in b else []) ) ) for k in set(a)|set(b) ] )
def sub(a,b) : return dict( filter(lambda (k,p):p, [(k, pairs( lumis(a[k]) - (lumis(b[k]) if k in b else set()))) for k in a ] ) )

jsons = dict((dir,dict( (f,eval(open(dir+'/'+f).readline())) for f in os.listdir(dir))) for dir in ['json','dsets_el','dsets_mu'])
for key,val in jsons.items() :
    print key
    for j1,j2 in itertools.combinations(val.keys(),2) :
        if intersect(val[j1],val[j2]) : print '\t', j1, j2, "intersect"
print
print '-'*20
print

for dir1,dir2 in itertools.combinations(jsons,2) :
    overlapping_pairs = [(j1,j2) for j1 in jsons[dir1] for j2 in jsons[dir2] if intersect(jsons[dir1][j1],jsons[dir2][j2])]
    print dir1,dir2
    print '\n'.join(str(p) for p in overlapping_pairs)
    print
print '-'*20
print


for dir1,dir2 in itertools.combinations(jsons,2) :
    if 'json'!=dir1 : continue
    tot1 = reduce(union, jsons[dir1].values() )
    tot2 = reduce(union, jsons[dir2].values() )

    print '%s - %s :'%(dir1,dir2)
    print sub(tot1,tot2)

    print '%s - %s :'%(dir2,dir1)
    print sub(tot2,tot1)
