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

    def __init__(self, calculateLumi = True, verbose = False) :
        self.lumisByRun = collections.defaultdict(list)
        self.calculateLumi = calculateLumi
        self.verbose = verbose
        self.moreName="see below"

    def uponAcceptance(self,eventVars) :
        self.lumisByRun[eventVars["run"]].append(eventVars["lumiSection"])
    
    def varsToPickle(self) : return ["lumisByRun"]

    def outputSuffix(self) : return ".json"

    def lumi(self, json) :
        return (-1 if not self.calculateLumi else
                utils.luminosity.recordedInvMicrobarns(json,False)/1e6 )

    def mergeFunc(self, products) :
        d = collections.defaultdict(list)
        for lumisByRun in products["lumisByRun"] :
            for run,lumis in lumisByRun.iteritems() :
                d[run] += lumis

        d2 = {}
        for run,lumis in d.iteritems() :
            d2[run] = sorted(set(lumis))
            for ls in d2[run] :
                if 1 < lumis.count(ls) and self.verbose :
                    print "Run %d ls %d appears %d times in the lumiTree."%(run,ls,lumis.count(ls))

        json = utils.jsonFromRunDict(d2)
        with open(self.outputFileName,"w") as file: print >> file, str(json).replace("'",'"')

        print "Wrote %.4f/pb json to : %s"%(self.lumi(json),self.outputFileName)
        print utils.hyphens
#####################################
class duplicateFileRunLumiCheck(analysisStep) :
    def __init__(self) :
        self.runLumis = collections.defaultdict(set)

    def uponAcceptance(self,ev) :
        self.runLumis[ev['treeFileName']].add((ev['run'],ev['lumiSection']))

    def varsToPickle(self) : return ["runLumis"]

    def mergeFunc(self, products) :
        for runLumis in products['runLumis'] : self.runLumis.update(runLumis)

        runLumis = sorted(sum([list(runLumis) for runLumis in self.runLumis.values()],[]))
        for key in sorted(set(runLumis)) :
            runLumis.remove(key)
        dups = set(runLumis)
        print "N duplicate lumis", len(dups)
        for dup in sorted(dups) :
            print dup, [f for f,runLumis in self.runLumis.iteritems() if dup in runLumis]
        print
        print sorted(dups)

class duplicateEventCheck(analysisStep) :
    def __init__(self, onlyKeys = []) :
        self.events = collections.defaultdict(list)
        self.onlyKeys = onlyKeys

    def uponAcceptance(self,ev) :
        key = (ev['run'],ev['lumiSection'])
        if key in self.onlyKeys or not self.onlyKeys:  self.events[key].append(ev["event"])

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
