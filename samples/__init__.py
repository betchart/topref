for module in [
    "mc",
    "mcOld",
    "jetmet",
    "ht",
    "muon",
    "electron",
    "photon",
    "signalSkim",
    "wpol",
    "mumu",
    "top",
    "ewk",
    "qcd",
    #
    "electron16",
    "muon16",
    "qcd16",
    "ewk16",
    "top16",
    #
    "photon17",
    "top17",
    ] : exec("from __%s__ import %s"%(module,module))
