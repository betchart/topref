for module in [
    "electron16",
    "muon16",
    "qcd112",
    "ewk112",
    "top112",
    ] : exec("from __%s__ import %s"%(module,module))
