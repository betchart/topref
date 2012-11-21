for module in [
    "electron16",
    "muon16",
    "qcd16",
    "ewk112",
    "top112",
    ] : exec("from __%s__ import %s"%(module,module))
