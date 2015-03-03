import supy
ewk119_fnal = supy.samples.SampleHolder()

eos = 'utils.fileListFromDisk(location="/eos/uscms/store/user/bbetchar/'


#dynj = "2013_01_25_01_48_16/DY%dJetsToLL_M-50_TuneZ2Star_8TeV-madgraph.Summer12_DR53X-PU_S10_START53_V7%s-v1.AODSIM"
dynj = "v119/DY%dJetsToLL_M-50_TuneZ2Star_8TeV-madgraph.Summer12_DR53X-PU_S10_START53_V7%s-v1.AODSIM"
ewk119_fnal.add("dy1j_mg", '%s/%s")'%(eos,dynj%(1,'A')), xs = 561.0, nCheck=24045248)
ewk119_fnal.add("dy2j_mg", '%s/%s")'%(eos,dynj%(2,'C')), xs = 181.0, nCheck=21852156)
ewk119_fnal.add("dy3j_mg", '%s/%s")'%(eos,dynj%(3,'A')), xs =  51.1, nCheck=11015445)
ewk119_fnal.add("dy4j_mg", '%s/%s")'%(eos,dynj%(4,'A')), xs =  23.04, nCheck=6402827)


#wnj = "2013_01_25_01_32_38/W%dJetsToLNu_TuneZ2Star_8TeV-madgraph.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM"
wnj = "v119/W%dJetsToLNu_TuneZ2Star_8TeV-madgraph.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM"
ewk119_fnal.add("w1j_mg", '%s/%s")'%(eos,wnj%1), xs = 5400.0 , nCheck=23141598)
ewk119_fnal.add("w2j_mg", '%s/%s")'%(eos,wnj%2), xs = 1750.0 , nCheck=34044921)
ewk119_fnal.add("w3j_mg", '%s/%s")'%(eos,wnj%3), xs =  519.0 , nCheck=15539503)
ewk119_fnal.add("w4j_mg", '%s/%s")'%(eos,wnj%4), xs =  214.0 , nCheck=13382803)

ewk119_fnal.add("wbb_mg", '%s/%s")'%(eos,'TOP/automated/2013_11_25_22_33_48'), xs = 211.3 , nCheck=16267456)
