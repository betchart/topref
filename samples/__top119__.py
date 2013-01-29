from supy.samples import SampleHolder
from supy.sites import srm
top119 = SampleHolder()

srm_burt = srm + '/bbetchar/TOP/automated'
singleT = srm_burt + '/2013_01_25_01_59_33/%s_%s_TuneZ2star_8TeV-powheg-tauola.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM/")'
ttj_ph = '2013_01_25_01_25_48/TT_CT10_TuneZ2star_8TeV-powheg-tauola.Summer12_DR53X-PU_S10_START53_V7A-v%d.AODSIM'

top119.add("ttj_ph", '+'.join(['%s/%s")'%(srm_burt,ttj_ph%v) for v in [1,2]]), xs = 211.0 )
#top119.add("ttj_mn", '%s/2013_01_15_05_49_45/")'%srm_burt, xs = 211.1 )

top119.add("top_s_ph", singleT%('T','s-channel'), xs = 2.82 )
top119.add("top_t_ph", singleT%('T','t-channel'), xs = 47.0 )
top119.add("top_tW_ph", singleT%('T','tW-channel-DR'), xs = 10.7 )
top119.add("tbar_s_ph", singleT%('Tbar','s-channel'), xs = 1.57 )
top119.add("tbar_t_ph", singleT%('Tbar','t-channel'), xs = 25.0 )
top119.add("tbar_tW_ph", singleT%('Tbar','tW-channel-DR'), xs = 10.7 )

top119.add("syncMC", '''eval('["/vols/cms04/bbetchar/sync/mc53X_v119.root"]')''',       xs = 211 )
top119.add("syncMC_db", '''eval('["/vols/cms04/bbetchar/sync/mc53X_v119_dB.root"]')''', xs = 211 )
