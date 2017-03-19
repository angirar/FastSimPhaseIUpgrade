# To install
cmsrel CMSSW_9_0_0_pre6
cd CMSSW_9_0_0_pre6/src
cmsenv
git clone https://github.com/angirar/FastSimPhaseIUpgrade.git FastSimulation
scram b -rj32

# To run
create a file with generated events-
source FastSimulation/SimplifiedGeometryPropagator/test/gen.sh
pass the generated events to simulation-
cmsRun FastSimulation/SimplifiedGeometryPropagator/python/conf_cfg.py
# to run validation do instead
cmsRun FastSimulation/SimplifiedGeometryPropagator/python/conf_validation_cfg.py