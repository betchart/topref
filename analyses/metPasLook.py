#!/usr/bin/env python

import analysis,utils
        
a=analysis.analysis(name="metPasLook",
                    outputDir="/vols/cms02/elaird1/tmp/"
                    )

a.addSampleSpec(sampleName="JetMETTau.Run2010A-Jul6thReReco_v1.RECO_cleanEvent",
                listOfFileNames=utils.fileListFromDisk(location="/vols/cms02/elaird1/05_skims/JetMETTau.Run2010A-Jul6thReReco_v1.RECO",pruneList=False,nMaxFiles=-1),
                listName="metPasLook",
                isMc=False,
                nEvents=-1)

#a.addSampleSpec(sampleName="JetMETTau.Run2010A-Jun14thReReco_v2.RECO.Bryn",
#                listOfFileNames=utils.fileListFromSrmLs(location="/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/bm409//ICF/automated/2010_07_20_16_52_06/",nMaxFiles=4),
#                listName="metPasLook",
#                isMc=False,
#                nEvents=-1)

a.go(loop=True,
     plot=True,
     profile=False,
     nCores=6,
     splitJobsByInputFile=True
     )
