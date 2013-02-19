This project is a satellite of [supy](github.com/elaird/supy) for [TopRef](http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/Betchart/TopRefTuple/) format TTrees.

# License
GPLv3 (http://www.gnu.org/licenses/gpl.html)

# Dependencies
[supy](github.com/elaird/supy) dependencies: [ROOT](http://root.cern.ch) (>=5.27.06) and python (2.x, x>=6).
CMSSW is not required, but is often the easiest way to get functioning PyROOT.

# Quickstart
1. Set up ROOT and Python.  For example, with CMSSW
```bash
cd /<somepath>/CMSSW_5_3_8/src && cmsenv && cd -
```
2. Initialize local repository (only once)
```bash
git clone git://github.com/betchart/topref.git   # no repo write permission, or
#git clone git://github.com/<username>/topref.git # if you have forked it
cd topref/
git submodule update --init                      # Initialize the supy submodule
```
3. Configure path variables
```bash
export PYTHONPATH=$PYTHONPATH:`pwd`              # Add supy/ parent directory to python path
export PATH=$PATH:`pwd`/supy/bin                 # optionally add supy/bin to your path
```
4. Run Tests
```bash
supy-test
```

# Bugs
Please report problems at [issues](https://github.com/betchart/topref/issues)
