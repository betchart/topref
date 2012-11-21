from supy.samples import SampleHolder
from supy.sites import srm
qcd112 = SampleHolder()

srm_burt = srm+'/bbetchar/TOP/automated'

em = "2012_11_19_06_13_24/QCD_Pt_%s_EMEnriched_TuneZ2star_8TeV_pythia6.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM"
em_xs = [
    ( 20, 2.886e8),
    ( 30, 7.433e7),
    ( 80, 1191000.0),
    (170,   30990.0),
    (250,    4250.0),
    (350,     810.0),
    (None, None)
    ]
for (lo,xs),(hi,_) in zip(em_xs[:-1],em_xs[1:]) :
    bin = "_".join(str(i) for i in  [lo,hi] if i!=None)
    qcd112.add("qcd_em_%s"%bin, '%s/%s")'%(srm_burt,em%bin), xs)

mu = "2012_11_18_06_05_40/QCD_Pt-%s_MuEnrichedPt5_TuneZ2star_8TeV_pythia6.Summer12_DR53X-PU_S10_START53_V7A-v%d.AODSIM"
mu_xs = [
    (  15, 7.022e8,    2),
    (  20, 2.870e8,    1),
    (  30, 6.609e7,    1),
    (  50, 8082000.0,  1),
    (  80, 1024000.0,  1),
    ( 120,  157800.0,  1),
    ( 170,   34020.0,  1),
    ( 300,    1757.0,  1),
    ( 470,     115.2,  1),
    ( 600,      27.01, 1),
    ( 800,       3.57, 1),
    (1000,       0.774,1),
    (None, None, None),
    ]
for (lo,xs,v),(hi,_,_) in zip(mu_xs[:-1],mu_xs[1:]) :
    bin = "to".join(str(i) for i in  [lo,hi] if i!=None)
    qcd112.add("qcd_mu_%s"%(bin.replace("to","_")), '%s/%s")'%(srm_burt,mu%(bin,v)), xs)
