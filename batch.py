import sys,os
from optparse import OptionParser

parser = OptionParser("usage: %prog tagSuffix [options]")
parser.add_option("--all",   dest = "doAll",     default = False,  action='store_true',  help = "submit all samples")
parser.add_option("--tops",  dest = "justTops",  default = False,  action='store_true',  help = "submit ttbar samples")
parser.add_option("--wjets", dest = "justWJets", default = False,  action='store_true',  help = "submit wjets samples")
parser.add_option("--nodata",dest = "noData",    default = False,  action='store_true',  help = 'submit all samples except data')

options,args = parser.parse_args()
if len(args)!=1 or sum(int(i) for i in [options.doAll,options.justTops,options.justWJets,options.noData])!=1:
    parser.print_help()
    exit()

tagSuffix=args[0]
tags = tagSuffix if ('top_' in tagSuffix or 'QCD_' in tagSuffix) else ','.join(pre+tagSuffix for pre in ['top_','QCD_'])
lep = 'El' if '_el_' in tags else 'Mu' if '_mu_' in tags else ''
top = 'ph' if '_ph_' in tags else 'mn' if '_mn_' in tags else ''
assert lep

groupsAll = {80:[ '%s.B.1.jw'%lep,
                  '%s.C.2.jw'%lep,
                  '%s.D.1.jw'%lep,
                  'w2j_mg.pu',
                  'ttj_%s.wGG.pu'%top,
                  ],
             10:['top_s_ph.pu',
                 'top_t_ph.pu',
                 'top_tW_ph.pu',
                 'tbar_s_ph.pu',
                 'tbar_t_ph.pu',
                 'tbar_tW_ph.pu']
             }

groupsTop = {142:['ttj_%s.wGG.pu'%top,],
             80:['ttj_%s.wQG.pu'%top,'ttj_%s.wQQ.pu'%top,],
             30:['ttj_%s.wAG.pu'%top,]
             }

groupsW = {80:['w2j_mg.pu'],
          40:['w%dj_mg.pu'%d for d in [1,3,4]]
          }


groupsNoData = {20:['.'.join([lep,part,'jw']) for part in ['A.1','A.2','B.1','C.1','C.2','C.3','D.1'] for lep in ['Mu','El']]}

cmds = []
groups = groupsW if options.justWJets else groupsTop if options.justTops else groupsNoData if options.noData else groupsAll
samplesOption = '--omit' if options.noData else '--samples'
for nslices,samples in groups.items() :
    cmds.append('supy analyses/topAsymm.py --tags %s %s %s --loop 2 --slices %d --batch'%(tags,samplesOption,','.join(samples),nslices))
if options.doAll : cmds.append('supy analyses/topAsymm.py --tags %s --omit %s --loop 2 --slices %d --batch'%(tags,','.join(sum([samples for samples in groups.values()],[])),40))

for cmd in cmds : print cmd

query = raw_input('Do it?')
if query in [1,True,'y','yes','Y','Yes','YES'] :
    for cmd in cmds :
        os.system(cmd)