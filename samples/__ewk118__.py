from supy.samples import SampleHolder
from supy.sites import srm
ewk118 = SampleHolder()

srm_burt = srm+'/bbetchar/TOP/automated'

dynj = "2013_01_15_05_33_36/DY%dJetsToLL_M-50_TuneZ2Star_8TeV-madgraph.Summer12_DR53X-PU_S10_START53_V7%s-v1.AODSIM"
ewk118.add("dy1j_mg", '%s/%s", alwaysUseLastAttempt = True)'%(srm_burt,dynj%(1,'A')), xs = 561.0)
ewk118.add("dy2j_mg", '%s/%s")'%(srm_burt,dynj%(2,'C')), xs = 181.0)
ewk118.add("dy3j_mg", '%s/%s")'%(srm_burt,dynj%(3,'A')), xs =  51.1)
ewk118.add("dy4j_mg", '%s/%s")'%(srm_burt,dynj%(4,'A')), xs =  23.04)


wnj = "2013_01_15_05_20_37/W%dJetsToLNu_TuneZ2Star_8TeV-madgraph.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM"
ewk118.add("w1j_mg", '%s/%s")'%(srm_burt,wnj%1), xs = 5400.0 )
ewk118.add("w2j_mg", '%s/%s")'%(srm_burt,wnj%2), xs = 1750.0 )
ewk118.add("w3j_mg", '%s/%s")'%(srm_burt,wnj%3), xs =  519.0 )
ewk118.add("w4j_mg", '%s/%s")'%(srm_burt,wnj%4), xs =  214.0 )
