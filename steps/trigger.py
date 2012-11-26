import math,collections,re, ROOT as r
from supy import analysisStep,utils
#####################################
class hltKeys(analysisStep) :
    def setup(self,*_) : self.keys = collections.defaultdict(lambda : collections.defaultdict(set) )
    def uponAcceptance(self, ev) :
        d = self.keys[str(ev['hltKey'])]
        for t,p in ev['prescaled'] :
            if not any(item in t for item in ['Central','HLT_Ele25_CaloIdV']) : continue
            else: d[str(t)].add( p )
    def varsToPickle(self) : return ['keys']
    def endFunc(self,*_) : self.keys = dict(self.keys)
    def mergeFunc(self, products) :
        allkeys = collections.defaultdict(lambda : collections.defaultdict(set) )
        for keys in products['keys'] :
            for key,ps in keys.items() :
                for t,p in ps.items() :
                    allkeys[key][t] |= p
        with open(self.outputFileName,"w") as outFile : 
            for key in sorted(allkeys) :
                print>>outFile, key
                for t,ps in allkeys[key].items() :
                    if 1 in ps :
                        print>>outFile, '\t',len(ps), t, ps
        print "Wrote ", self.outputFileName
    def outputSuffix(self) :
        return "_hltkeys.txt"
#####################################
class NameDump(analysisStep) :

    def __init__(self,triggerLevel = ""):
        self.varName = triggerLevel + "triggered"
        self.moreName = self.varName

    def select (self,eventVars) :
        for pair in eventVars[self.varName] :
            print pair.first
        return True
#####################################
class Counts(analysisStep) :

    def __init__(self, useCache = False) :
        self.useCache = useCache
        self.counts = collections.defaultdict(int)
        self.cached = {}

    def uponAcceptance(self, eventVars) :
        if not self.useCache :
            for pair in eventVars["triggered"] :
                if pair.second : self.counts[pair.first] += 1
        else :
            key = (eventVars["run"], eventVars["lumiSection"])
            m = eventVars["triggered"]
            self.updateCached(key, m)
            for t in self.cached[key] :
                if m[t] : self.counts[t] += 1

    def updateCached(self, key, m) :
        if key not in self.cached :
            self.cached[key] = [pair.first for pair in m]

    def outputSuffix(self) :
        return "_triggerCounts.txt"

    def varsToPickle(self) :
        return ["counts"]

    def mergeFunc(self, products) :
        def mergedCounts(l) :
            out = collections.defaultdict(int)
            for d in l :
                for key,value in d.iteritems() :
                    out[key] += value
            return out

        counts = mergedCounts(products["counts"])
        outFile = open(self.outputFileName,"w")

        maxNameLength = max([len(key) for key in counts.keys()])
        maxCountLength = max([len(str(value)) for value in counts.values()])
        for key in sorted(counts.keys()) :
            outFile.write("%s    %s\n"%(key.ljust(maxNameLength), str(counts[key]).ljust(maxCountLength)))
        outFile.close()
        print "The trigger counts file %s has been written."%self.outputFileName
        print utils.hyphens
#####################################
class lowestUnPrescaledTriggerFilter(analysisStep) :
    def select (self,ev) : return ev["lowestUnPrescaledTrigger"] is not None
#####################################
class lowestUnPrescaledTriggerHistogrammer(analysisStep) :
    def __init__(self, collectVersions = True, drop = []) :
        self.key = "lowestUnPrescaledTrigger"
        self.collectVersions = collectVersions
        self.drop = drop

    def makeLabels(self, eventVars) :
        self.label = {}
        paths = dict.__getitem__(eventVars,self.key).sortedListOfPaths
        self.labels = []
        for path in paths :
            s = path.split("_")
            if self.collectVersions and s[0]=="HLT" and s[-1][0]=="v" :
                self.label[path] = "_".join(s[1:-1])
            else :
                self.label[path] = path
            for item in self.drop : self.label[path] = self.label[path].replace(item,'.')
            if not self.labels or self.labels[-1]!=self.label[path] :self.labels.append(self.label[path])
        self.nBins = len(self.labels)
        
    def uponAcceptance(self, eventVars) :
        if not hasattr(self, "labels") :
            self.makeLabels(eventVars)

        iBin = self.labels.index(self.label[eventVars[self.key]])
        self.book.fill(iBin, self.key, self.nBins, 0.0, self.nBins, title = ";lowest un-prescaled path;events / bin", xAxisLabels = self.labels)
#####################################
class hltFilter(analysisStep) :
    def __init__(self,hltPathName):
        self.hltPathName = hltPathName
        self.moreName = self.hltPathName

    def select (self,eventVars) :
        return eventVars["prescaled"][self.hltPathName]
#####################################
class hltFilterList(analysisStep) :

    def __init__(self,hltPathNames):
        self.hltPathNames = hltPathNames
        self.moreName = "any of "+str(self.hltPathNames)

    def select (self,ev) :
        return any( ev[path] for path in self.hltPathNames )
#####################################
class hltPrescaleHistogrammer(analysisStep) :

    def __init__(self,listOfHltPaths) :
        self.listOfHltPaths = listOfHltPaths
        self.moreName = ','.join(self.listOfHltPaths).replace("HLT_","")
        self.nBinsX = len(self.listOfHltPaths)
        self.key = "HltPrescaleHisto"

    def uponAcceptance(self,eventVars) :
        for iPath in range(len(self.listOfHltPaths)) :
            value = eventVars["prescaled"][self.listOfHltPaths[iPath]]
            if value<=0.0 : continue
            self.book.fill( (iPath,math.log10(value)), self.key, (self.nBinsX,100), (-0.5,-0.5), (self.nBinsX-0.5,4,5),
                            title="hltPrescaleHisto;;log_{10}(prescale value);events / bin", xAxisLabels = self.listOfHltPaths)
#####################################
class hltTurnOnHistogrammer(analysisStep) :

    def __init__(self, var = None, binsMinMax = None, probe = None, tags = None, permissive = False) :
        for item in ["var","binsMinMax","probe","tags","permissive"] :
            setattr(self,item,eval(item))

        tags = "{%s}"%(','.join([t.replace("HLT_","") for t in self.tags]))
        probe = self.probe.replace("HLT_","")
        
        self.tagTitle   = ( "tag_%s_%s_%s" % (probe,tags,var), "pass %s;%s;events / bin" % (tags,var))
        self.probeTitle = ( "probe_%s_%s_%s" % (probe,tags,var), "pass %s given %s;%s;events / bin" % (probe, tags, var) )
        self.effTitle   = ( "turnon_%s_%s_%s" % (probe,tags,var), "%s Turn On;%s;P(%s | %s)" % (probe,var,probe,tags))

        self.moreName = "%s by %s, given %s;" % (probe, var, tags)

    def uponAcceptance(self,eventVars) :
        if not any([eventVars["prescaled"][t] for t in self.tags]) : return
        if (not self.permissive) and 1 != eventVars["prescaled"][self.probe] : return
        
        for t in ([self.tagTitle] if not eventVars["prescaled"][self.probe] else [self.tagTitle,self.probeTitle]) :
            self.book.fill( eventVars[self.var], t[0], self.binsMinMax[0],self.binsMinMax[1],self.binsMinMax[2], title = t[1] )
    
    def mergeFunc(self, products) :
        probe = r.gDirectory.Get(self.probeTitle[0])
        tag = r.gDirectory.Get(self.tagTitle[0])
        if not (tag and probe) : return
        efficiency = probe.Clone(self.effTitle[0])
        efficiency.SetTitle(self.effTitle[1])
        efficiency.Divide(probe,tag,1,1,"B")
        for bin in [0,self.binsMinMax[0]+1] :
            efficiency.SetBinContent(bin,0)
            efficiency.SetBinError(bin,0)
        efficiency.Write()
        print "Output updated with TurnOn %s."%self.effTitle[0]
#####################################
class triggerScan(analysisStep) :
    def __init__(self, pattern = r".*", prescaleRequirement = "True", tag = "") :
        self.tag = tag
        self.pattern = pattern
        self.prescaleRequirement = prescaleRequirement
        self.moreName = "%s; %s"%(pattern, prescaleRequirement)

        self.triggerNames = collections.defaultdict(set)
        self.counts = collections.defaultdict(int)

    def uponAcceptance(self,eventVars) :
        key = (eventVars["run"],eventVars["lumiSection"])
        self.counts[key] += 1
        if key not in self.triggerNames :
            for name,prescale in eventVars["prescaled"] :
                if re.match(self.pattern,name) and eval(self.prescaleRequirement) :
                    self.triggerNames[key].add(name)

    def varsToPickle(self) : return ["triggerNames","counts"]
        
    def mergeFunc(self,products) :
        def update(a,b) : a.update(b); return a;
        self.triggerNames = reduce(update, products["triggerNames"], dict())
        self.counts = reduce(update, products["counts"], dict())

        names = sorted(list(reduce(lambda a,b: a|b, self.triggerNames.values(), set())))

        reducedNames = []
        reducedCounts = []
        order = sorted(self.triggerNames.keys())
        for key in order :
            if (not reducedNames) or reducedNames[-1] != self.triggerNames[key] :
                reducedNames.append(self.triggerNames[key])
                reducedCounts.append(self.counts[key])
            else : reducedCounts[-1] += self.counts[key]

        hist = r.TH2D(self.tag if self.tag else "triggerScan",
                      "(%s) with (%s);epochs of %s;;events"%(self.pattern, self.prescaleRequirement, self._analysisStep__outputFileStem.split('/')[-1]),
                      len(reducedNames),0,len(reducedNames),len(names),0,len(names))
        for i,name in enumerate(names) : hist.GetYaxis().SetBinLabel(i+1,name.replace("HLT_",""))
        for i in range(len(reducedNames)) : hist.GetXaxis().SetBinLabel(i+1,"%d"%(i+1))
        
        for i, iNames, iCount in zip(range(len(reducedNames)), reducedNames, reducedCounts) :
            for name in iNames :
                j = names.index(name)
                hist.SetBinContent(i+1,j+1,iCount)

        hist.Write()
        print "Output updated with triggerScans %s."%self.tag
#####################################
class prescaleScan(analysisStep) :
    '''Currently broken due to double duty "prescaled" branch (topRef/tree).'''
    def __init__(self, trigger = None, ptMin = None, triggeringPt = "") :
        self.trigger = trigger
        self.triggeringPt = triggeringPt
        self.ptMin = ptMin
        self.moreName = "%s; %.1f<pt mu"%(trigger, ptMin)

    def uponAcceptance(self,ev) :
        
        if not self.ptMin < ev[self.triggeringPt] : return
        prescale = ev["prescaled"][self.trigger]
        if not prescale : return
        name = self.trigger+"_p%d"%prescale
        self.book.fill( ev['triggered'][self.trigger], name, 2,0,2, title = '%s;Fail / Pass;event / bin'%(name))
#####################################
class anyTrigger(analysisStep) :
    def __init__(self, sortedListOfPaths = [], unreliable = {}) :
        self.sortedListOfPaths = sortedListOfPaths
        self.unreliable = unreliable
        
        self.moreName = "any of "+utils.contract(self.sortedListOfPaths)
        
    def select(self, ev) :
        return any(ev['triggered'][item] for item in self.sortedListOfPaths if item not in self.unreliable or ev['prescaled'][item] not in self.unreliable[item])

#####################################
class prescaleLumiEpochs(analysisStep) :
    '''Currently broken due to double duty "prescaled" branch (topRef/tree).'''

    def __init__(self, triggers) :
        self.triggers = zip(*triggers)[0]
        self.thresh = zip(*triggers)[1]

    def setup(self,ignored1,ignored2) :
        dd = collections.defaultdict
        self.epochLumis = dd(lambda : dd(set)) #set of lumis by run by (prescale set)
        self.lumis = dd(set) #set of lumis by run
    
    def uponAcceptance(self,ev) :
        run,lumi = (ev["run"],ev["lumiSection"])
        if lumi in self.lumis[run] : return
        key = tuple([ev["prescaled"][trigger] for trigger in self.triggers])
        self.epochLumis[key][run].add(lumi)
        self.lumis[run].add(lumi)
        
    def varsToPickle(self) : return ["staticEpochLumis"]

    @property
    def staticEpochLumis(self) : return dict((epoch,dict(self.epochLumis[epoch].iteritems())) for epoch in self.epochLumis)

    def mergeFunc(self,products) :
        self.setup(None,None)
        for eLumis in products["staticEpochLumis"] :
            for epoch in eLumis :
                for run in eLumis[epoch] :
                    self.epochLumis[epoch][run] |= eLumis[epoch][run]

        lumiProbs = []
        with open(self.outputFileName,"w") as file :
            print >> file, self.triggers
            print >> file, ""
            for epoch in self.epochLumis:
                json = utils.jsonFromRunDict(self.epochLumis[epoch])
                lumi = utils.luminosity.recordedInvMicrobarns(json)/1e6
                probs = [(thresh,1./prescale) for thresh,prescale in zip(self.thresh,epoch)[:epoch.index(1.0)+1] if prescale]
                inclusive = [(probs[i][0],probs[i+1][0] if (i<len(probs)-1) else None,utils.unionProbability(zip(*probs[:i+1])[1])) for i in range(len(probs))]
                lumiProbs.append((lumi,inclusive))
                print >> file, epoch
                print >> file, json
                print >> file, lumi, " /pb"
                print >> file, [(thresh,"%.4f"%w) for thresh,nextThresh,w in inclusive]
                print >> file, ""
        print "Wrote prescale jsons by epoch : "+self.outputFileName

        lumiVal = sum(lumi for lumi,probs in lumiProbs)
        lumiHist = r.TH1D("susyTreeLumi","luminosity from susyTree;;1/pb",1,0,1)
        lumiHist.Fill(0.5,lumiVal)
        lumiHist.Write()
        
        ptMax = 10+max(self.thresh)
        weight = r.TH1D("prescaleWeight","prescaleWeight",ptMax,0,ptMax)
        for lumi,probs in lumiProbs :
            epochWeight = lumi / lumiVal
            for thresh,nextThresh,prob in probs:
                for bin in range(thresh,nextThresh if nextThresh else ptMax) :
                    weight.Fill(bin, epochWeight * prob)
        weight.Write()
