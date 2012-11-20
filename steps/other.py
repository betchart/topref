import array,collections, ROOT as r
from supy import analysisStep,utils,steps
#####################################
class handleChecker(analysisStep) :
    def __init__(self, matches = ["handle", "Handle", "valid", "Valid"]) :
        self.run = False
        self.matches = matches

    def ofInterest(self, eventVars) :
        out = []
        for var in eventVars :
            for item in self.matches :
                if item in var :
                    out.append(var)
        return list(set(out))
    
    def uponAcceptance(self, eventVars) :
        if self.run : return
        self.run = True

        true = []
        false = []
        for var in self.ofInterest(eventVars) :
            if eventVars[var] :
                true.append(var)
            else :
                false.append(var)
        print "True:",sorted(true)
        print
        print "False:",sorted(false)
#####################################
class jsonMaker(analysisStep) :

    def __init__(self, calculateLumi = True, pixelLumi = True, debug = False) :
        self.lumisByRun = collections.defaultdict(list)
        self.calculateLumi = calculateLumi
        self.pixelLumi = pixelLumi
        self.debug = debug
        self.moreName="see below"

    def uponAcceptance(self,eventVars) :
        self.lumisByRun[eventVars["run"]].append(eventVars["lumiSection"])
    
    def varsToPickle(self) : return ["lumisByRun"]

    def outputSuffix(self) : return ".json"

    def lumi(self, json) :
        if not self.calculateLumi : return -1.0
        if self.pixelLumi :
            return utils.luminosity.recordedInvMicrobarns(json)/1e6
        else :
            dct = utils.getCommandOutput("lumiCalc2.py overview -i %s"%self.outputFileName)
            assert not dct["returncode"],dct["returncode"]
            assert not dct["stderr"],dct["stderr"]
            s = dct["stdout"]
            if self.debug : print s[s.find("Total"):]
            m = "Recorded(/"
            i = s.rindex(m) + len(m)
            units = s[i-1:i+2]
            factor = {"/fb":1.0e3, "/pb":1.0, "/nb":1.0e-3, "/ub":1.0e-6}
            assert units in factor,units
            i2 = dct["stdout"].rindex("|")
            i1 = dct["stdout"][:i2].rindex("|")
            return float(dct["stdout"][1+i1:i2])*factor[units]

    def mergeFunc(self, products) :
        d = collections.defaultdict(list)
        for lumisByRun in products["lumisByRun"] :
            for run,lumis in lumisByRun.iteritems() :
                d[run] += lumis

        d2 = {}
        for run,lumis in d.iteritems() :
            d2[run] = sorted(set(lumis))
            for ls in d2[run] :
                if 1 < lumis.count(ls) :
                    print "Run %d ls %d appears %d times in the lumiTree."%(run,ls,lumis.count(ls))

        json = utils.jsonFromRunDict(d2)
        with open(self.outputFileName,"w") as file: print >> file, str(json).replace("'",'"')

        print "Wrote %.4f/pb json to : %s"%(self.lumi(json),self.outputFileName)
        print utils.hyphens
#####################################
class duplicateEventCheck(analysisStep) :
    def __init__(self) :
        self.events = collections.defaultdict(list)

    def uponAcceptance(self,ev) :
        self.events[(ev["run"], ev["lumiSection"])].append(ev["event"])

    def varsToPickle(self) : return ["events"]

    def mergeFunc(self, products) :
        def mergedEventDicts(l) :
            out = collections.defaultdict(list)
            for d in l :
                for key,value in d.iteritems() :
                    out[key] += value
            return out

        def duplicates(l) :
            s = set(l)
            for item in s :
                l.remove(item)
            return list(set(l))

        anyDups = False
        events = mergedEventDicts(products["events"])
        for runLs in sorted(events.keys()) :
            d = duplicates(events[runLs])
            if d :
                print "DUPLICATE EVENTS FOUND in run %d ls %d: %s"%(runLs[0], runLs[1], d)
                anyDups = True
        if not anyDups :
            print "No duplicate events were found."
#####################################






#####################################
class productGreaterFilter(analysisStep) :

    def __init__(self, threshold, variables, suffix = ""):
        self.threshold = threshold
        self.variables = variables
        self.moreName = "%s>=%.3f %s" % ("*".join(variables),threshold,suffix)

    def select (self,eventVars) :
        product = 1
        for var in self.variables : product *= eventVars[var]
        return product >= self.threshold
#####################################
class pickEventSpecMaker(analysisStep) :
    #https://twiki.cern.ch/twiki/bin/viewauth/CMS/WorkBookPickEvents

    def __init__(self) :
        self.events = []
        
    def outputSuffix(self) :
        return "_pickEvents.txt"
    
    def uponAcceptance(self, eventVars) :
        self.events.append( (eventVars["run"], eventVars["lumiSection"], eventVars["event"]) )
        
    def varsToPickle(self) :
        return ["events"]

    def mergeFunc(self, products) :
        out = open(self.outputFileName, "w")
        for events in products["events"] :
            for event in events :
                out.write("%14d:%6d:%14d\n"%event)
        out.close()
        print "The pick events spec. file %s has been written."%self.outputFileName

#####################################
class cutSorter(analysisStep) :
    def __init__(self, listOfSteps = [], applySelections = True ) :
        self.selectors = filter(lambda s: s.isSelector, listOfSteps)
        self.applySelections = applySelections
        self.moreName = "Applied" if applySelections else "Not Applied"
        self.bins = 1 << len(self.selectors)

    def select(self,eventVars) :
        selections = [s.select(eventVars) for s in self.selectors]
        self.book.fill( utils.intFromBits(selections), "cutSorterConfigurationCounts", 
                                  self.bins, -0.5, self.bins-0.5, title = ";cutConfiguration;events / bin")
        return (not self.applySelections) or all(selections)
        
    def endFunc(self, chains) :
        bins = len(self.selectors)
        self.book.fill(1, "cutSorterNames", bins, 0, bins, title = ";cutName", xAxisLabels = [sel.__class__.__name__ for sel in self.selectors])
        self.book.fill(1, "cutSorterMoreNames", bins, 0, bins, title = ";cutMoreName", xAxisLabels = [sel.moreName for sel in self.selectors])
#####################################
