for module in [
    "lepton112",
    "qcd112",
    "ewk112",
    "top112",
    "top118",
    ] : exec("from __%s__ import %s"%(module,module))
