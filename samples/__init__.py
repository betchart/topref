for module in [
    "electron16",
    "muon16",
    "qcd16",
    "ewk16",
    "top16",
    ] : exec("from __%s__ import %s"%(module,module))
