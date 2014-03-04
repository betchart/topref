import supy
top119_fnal = supy.samples.SampleHolder()

pnfs_cre = supy.sites.pnfs().replace('lpcsusyra1','cerba17')
pnfs_bb = supy.sites.pnfs().replace('lpcsusyra1','bbetchar')

singleT = pnfs_cre + 'v119/%s_%s_TuneZ2star_8TeV-powheg-tauola.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM/")'
ttj_ph = 'v119/TT_CT10_TuneZ2star_8TeV-powheg-tauola.Summer12_DR53X-PU_S10_START53_V7A-v%d.AODSIM'

top119_fnal.add("ttj_ph", '+'.join(['%s/%s")'%(pnfs_cre,ttj_ph%v) for v in [1,2]]), xs = 211.0 )
top119_fnal.add("extra_ph", pnfs_bb + '/TOP/automated/2013_11_19_19_52_22")', xs = 9999 )

top119_fnal.add("top_s_ph", singleT%('T','s-channel'), xs = 2.82 )
top119_fnal.add("top_t_ph", singleT%('T','t-channel'), xs = 47.0 )
top119_fnal.add("top_tW_ph", singleT%('T','tW-channel-DR'), xs = 10.7 )
top119_fnal.add("tbar_s_ph", singleT%('Tbar','s-channel'), xs = 1.57 )
top119_fnal.add("tbar_t_ph", singleT%('Tbar','t-channel'), xs = 25.0 )
top119_fnal.add("tbar_tW_ph", singleT%('Tbar','tW-channel-DR'), xs = 10.7 )


top119_fnal.add("calib_mn", pnfs_bb + '/TOP/automated/2013_11_25_22_52_14/")', xs = 211.1)
ttj_mg = '/TOP/automated/2013_11_25_23_05_40/TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola.Summer12_DR53X-PU_S10_START53_V19_ext%d-v1.AODSIM'
top119_fnal.add("calib_mg", '+'.join(['%s/%s")'%(pnfs_bb,ttj_mg%v) for v in [1,2]]), xs = 53.2)

axigluons = "/TOP/automated/2014_01_30_00_49_16/"

L200 = "HeavyGluonToTT_left_M-200_TuneZ2star_8TeV-madgraph-tauola.Summer12_DR53X-PU_S10_START53_V19-v1.AODSIM"
A200 = "HeavyGluonToTT_axial_M-200_TuneZ2star_8TeV-madgraph-tauola.Summer12_DR53X-PU_S10_START53_V19-v1.AODSIM"
R200 = "HeavyGluonToTT_right_M-200_TuneZ2star_8TeV-madgraph-tauola.Summer12_DR53X-PU_S10_START53_V19-v1.AODSIM"
A2K = "HeavyGluonToTT_axial_M-2000_TuneZ2star_8TeV-madgraph-tauola.Summer12_DR53X-PU_S10_START53_V19-v1.AODSIM"
R2K = "HeavyGluonToTT_right_M-2000_TuneZ2star_8TeV-madgraph-tauola.Summer12_DR53X-PU_S10_START53_V19-v1.AODSIM"
ZP = "ZprimeTtoTTU_M-220_TuneZ2star_8TeV-madgraph-tauola.Summer12_DR53X-PU_S10_START53_V19-v1.AODSIM"

top119_fnal.add("calib_L200", pnfs_bb + axigluons + L200 + '")', xs = 211)
top119_fnal.add("calib_A200", pnfs_bb + axigluons + A200 + '")', xs = 211)
top119_fnal.add("calib_R200", pnfs_bb + axigluons + R200 + '")', xs = 211)
top119_fnal.add("calib_A2K", pnfs_bb + axigluons + A2K + '")', xs = 211)
top119_fnal.add("calib_R2K", pnfs_bb + axigluons + R2K + '")', xs = 211)
top119_fnal.add("calib_ZP", pnfs_bb + axigluons + ZP + '")', xs = 211)

