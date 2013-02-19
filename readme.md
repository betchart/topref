This project is a satellite of [supy](https://github.com/elaird/supy) for [TopRef](http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/Betchart/TopRefTuple/) format TTrees.

## License
[GPLv3](http://www.gnu.org/licenses/gpl.html)

## Dependencies
[supy](https://github.com/elaird/supy) dependencies: 
* Python (2.x, x>=6)
* [ROOT](http://root.cern.ch) (>=5.27.06) with PyROOT enabled

Further dependencies:
* NumPy and [SciPy](http://www.scipy.org)

CMSSW is not required, but is often the easiest way to get functioning PyROOT.

## Quickstart
Steps 1 and 3 need to be run any time you start in a new shell; step 2 only needs to be done once per machine.

1. Set up ROOT and Python, for example, by running `cmsenv` from an appropriate path.

2. Initialize local repository (only once)
```bash
    git clone git://github.com/betchart/topref.git    # no repo write permission, or
    #git clone git://github.com/<username>/topref.git # if you have forked it
    cd topref/
    git submodule update --init                       # initialize the supy submodule
```
3. Configure path variables
```bash
    export PYTHONPATH=$PYTHONPATH:`pwd`              # add parent directory of supy to python path
    export PATH=$PATH:`pwd`/supy/bin                 # optionally add supy/bin to your path
```
4. Run Tests
```bash
    supy-test
```

## Bugs
Please report problems at [issues](https://github.com/betchart/topref/issues)
