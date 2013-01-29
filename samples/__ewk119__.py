from supy.samples import SampleHolder
from supy.sites import srm
ewk119 = SampleHolder()

srm_burt = srm+'/bbetchar/TOP/automated'

dynj = "2013_01_25_01_48_16/DY%dJetsToLL_M-50_TuneZ2Star_8TeV-madgraph.Summer12_DR53X-PU_S10_START53_V7%s-v1.AODSIM"
ewk119.add("dy1j_mg", '%s/%s", alwaysUseLastAttempt = True)'%(srm_burt,dynj%(1,'A')), xs = 561.0)
ewk119.add("dy2j_mg", '%s/%s")'%(srm_burt,dynj%(2,'C')), xs = 181.0)
ewk119.add("dy3j_mg", '%s/%s")'%(srm_burt,dynj%(3,'A')), xs =  51.1)
ewk119.add("dy4j_mg", '%s/%s")'%(srm_burt,dynj%(4,'A')), xs =  23.04)


wnj = "2013_01_25_01_32_38/W%dJetsToLNu_TuneZ2Star_8TeV-madgraph.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM"
ewk119.add("w1j_mg", '%s/%s")'%(srm_burt,wnj%1), xs = 5400.0 )
ewk119.add("w2j_mg", '%s/%s")'%(srm_burt,wnj%2), xs = 1750.0 )
ewk119.add("w3j_mg", '%s/%s")'%(srm_burt,wnj%3), xs =  519.0 )
ewk119.add("w4j_mg", '%s/%s")'%(srm_burt,wnj%4), xs =  214.0 )
