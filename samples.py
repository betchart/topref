import collections, configuration, calculables, copy

def specify(names = [], overrideLumi = None, xsPostWeights = None, effectiveLumi = None, nFilesMax = None, nEventsMax = None, weights = [], color = 1, markerStyle = 1 ) :
    assert not (overrideLumi and type(names)==list)
    if type(names) != list : names = [names]
    if type(weights) != list : weights = [weights]
    samplespec = collections.namedtuple("samplespec", "name overrideLumi xsPostWeights effectiveLumi nFilesMax nEventsMax weights color markerStyle")
    return [samplespec(name,overrideLumi,xsPostWeights,effectiveLumi,nFilesMax,nEventsMax,weights,color,markerStyle) for name in names]
    
class SampleHolder(dict) :
    sample = collections.namedtuple("sample", "filesCommand xs lumi ptHatMin")
    def __init__(self) : self.inclusiveGroups = []

    @property
    def inclusiveNames(self) : return set(sum(self.inclusiveGroups,()))
    
    def update(self, other) :
        assert type(other) is type(self), "%s is not a SampleHolder" % str(type(other))
        for key in other : assert key not in self, "%s already specified" % key
        dict.update(self,other)
        for group in other.inclusiveGroups : self.addInclusiveGroup(group)

    def add(self, name, filesCommand = None, xs = None, lumi = None, ptHatMin = None) :
        assert lumi or xs,                      "Underspecified sample: %s"%name
        assert not (lumi and (xs or ptHatMin)), "Overspecified sample: %s"%name
        self[name] = self.sample(filesCommand, xs, lumi, ptHatMin)

    def addInclusiveGroup( self, samplesGroup) :
        for sample in samplesGroup :
            assert 1==samplesGroup.count(sample), "Duplicate sample %s in inclusive group."%sample
            assert sample in self, "Unknown sample %s"%s
            assert self[sample].ptHatMin, "ptHatMin unspecified for sample: %s"%s

        inter = self.inclusiveNames.intersection(set(samplesGroup))
        assert not inter, "Samples { %s } already in other group."%",".join(inter)
        self.inclusiveGroups.append( tuple(sorted(samplesGroup, key = lambda name: self[name].ptHatMin) ) )

    def manageInclusive(self, sampleSpecs = [], applyPostWeightXS = True) :
        inclusiveSpecs = filter(lambda ss: ss.name in self.inclusiveNames, sampleSpecs )
        for spec in inclusiveSpecs :
            group = filter(lambda g: spec.name in g, self.inclusiveGroups)[0]
            greaterPtHat = group[ 1 + group.index(spec.name) : ]

            if not greaterPtHat : continue
            if applyPostWeightXS :
                if spec.xsPostWeights : print "Warning: %s has already specified xsPostWeights : Unsafe to manage inclusive samples!"%spec.name
                if spec.weights : print "Warning: %s already has weights { %s } : Unsafe to manage inclusive samples if other weights can be None!"%(spec.name,','.join(["%s %s"%(w.name(),w.moreName) for w in spec.weights]))
            
            nextPtHat = self[greaterPtHat[0]]
            modArgs = dict([(item,getattr(spec,item)) for item in spec._fields])
            modArgs["names"] = modArgs["name"]
            del modArgs["name"]
            modArgs["weights"] = copy.deepcopy(modArgs["weights"]) + [calculables.Other.pthatLess( nextPtHat.ptHatMin ) ]
            if applyPostWeightXS : modArgs["xsPostWeights"] = self[spec.name].xs - nextPtHat.xs
            sampleSpecs[sampleSpecs.index(spec)] = samples.specify( **modArgs )[0]
            
        return sampleSpecs

        
for module in configuration.samplesFiles() :
    exec("from samples%s import *"%module)
