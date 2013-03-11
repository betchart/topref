import os,collections
import ROOT as r

dir = "/vols/cms04/bbetchar/tmp//topAsymm/"
tags = [item for item in os.listdir(dir) if any(i in item for i in ['_mu_','_el_']) and '.' not in item]

lengths = collections.defaultdict(dict)

for tag in tags :
    suffix = '_plots.root'
    plotfiles = [item.replace(suffix,'') for item in os.listdir(dir+'/'+tag) if suffix in item]
    for pname in plotfiles :
        tfile = r.TFile.Open(dir+'/'+tag+'/'+pname+suffix)
        counts = tfile.Get('master/counts')
        key = (round(counts.GetBinContent(1),1), round(counts.GetBinContent(2),1))
        if key not in lengths[pname] : lengths[pname][key] = set()
        lengths[pname][key].add(tag)
        tfile.Close()

for key in sorted(lengths) :
    print
    print key,
    if len(lengths[key])>1 : print
    for key2,val2 in sorted(lengths[key].items(), key = lambda ((x1,x2),_):x1+x2) :
        print '\t', sum(key2), key2, ','.join(sorted(val2)) if len(lengths[key])>1 else ''
