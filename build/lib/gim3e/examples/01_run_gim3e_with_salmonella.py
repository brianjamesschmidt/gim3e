# A sample script to illustrate the use of GIM3E and check the 
# results your solver returns
#
# The tests were run with the 32-bit CPLEX 12.5.0 
# solver, with CoBRApy 0.2.x and Python 2.7 on MacOSX 10.8.4.
# There is a free academic CPLEX license available from IBM.
#
# If running cplex in Python/OSX, be sure to set 32 bit python at the shell
# command line since IBM hasn't released a 64-bit version for OSX.
# For macs:
# export VERSIONER_PYTHON_PREFER_32_BIT=yes
#
# Note GLPK does not work very well for GIM3E MILP problems
# 
# Note Gurobi 5.1 appears to work, but note integer tolerances
# are different than for CPLEX 12.5.
#

from gim3e.core import gim3e
import pickle
from cobra.manipulation import initialize_growth_medium
from copy import deepcopy
from types import *
from  cobra import __version__ as cobra_version


gim3e_dir = gim3e.__file__
n_remove_chars = len('/core/gim3e.py')
gim3e_dir = gim3e_dir[:(-1 * (n_remove_chars))]
data_dir = gim3e_dir + "/data/"

# read in the pre-processed omics data
fp = open(data_dir + "transcriptomics_dict.pickle", "rb")
transcriptomics_dict = pickle.load(fp)
fp.close()

fp = open(data_dir + "metabolomics_dict.pickle", "rb")
metabolomics_dict = pickle.load(fp)
fp.close()

# Make models not informed by omics first by initializing the media
import pickle
fp = open(data_dir + "salmonella_gem.pickle", "rb")
STM = pickle.load(fp)
fp.close()

if cobra_version == '0.2.0':
    STM_LB = deepcopy(STM)
else:
    STM_LB = STM.copy()
initialize_growth_medium(STM_LB, 'LB')

## This is how to set up LPM
# if cobra_version == '0.2.0':
#     STM_LPM = deepcopy(STM)
# else:
#     STM_LPM = STM.copy()
# initialize_growth_medium(STM_LPM, 'LPM')
## Protons should only have a net positive diffusion 
## from the extracellular space into the perplasm in LPM condition
# STM_LPM.reactions.get_by_id('Htex').lower_bound = 0

# Required numeric tolerances were 
# determined previously from FVA 
# on model reactions
selected_tolerance = 1E-8
selected_growth = 0.99
selected_penalty = 1.01

# Note detected metabolites that paired
# with the model are included even if
# they are blocked.  The algorithm
# automatically detects and warns
# about blocked metabolites.
# Note we can also add in virtual metabolites
# to merely monitor metabolites by formulating
# the model and not imposing flux requirements
gim3e_model, gim3e_FVA, total_penalty = gim3e.gim3e(STM_LB, expression_dict = transcriptomics_dict["LB"],
    expression_threshold = max(transcriptomics_dict["LB"].values()),
    metabolite_list = metabolomics_dict["lb"], 
    fraction_growth = selected_growth,
    relative_penalty_bound = selected_penalty,
    solver_tolerance = selected_tolerance,
    metabolite_flux_requirement = True,
    monitor_all_cellular_metabolites = True,
    MILP_formulation = True,
    run_FVA = True, reduce_model = False,
    trim_model = False, solver = 'cplex')

# Check the results against a reference run.  
filename = "STM_LB_FVA_reference.pickle"
fp = open(data_dir + filename, "rb")
gim3e_FVA_ref = pickle.load(fp)
fp.close()

# This is an optional conversion, but don't need to do just to compare solver performance
new_FVA = gim3e_FVA
reference_FVA = gim3e_FVA_ref
delta_dict_max = {}
delta_dict_min = {}
unsolved_min_one = {}
unsolved_max_one = {}
x = []
y = []
mismatch_dict = {}
shared_keys = set(new_FVA.keys()).intersection(set(reference_FVA.keys()))
for the_key in shared_keys:
    if ((not isinstance(new_FVA[the_key]['maximum'],NoneType)) & (not isinstance(reference_FVA[the_key]['maximum'],NoneType))):
        delta_dict_max[the_key] = abs(new_FVA[the_key]['maximum'] - reference_FVA[the_key]['maximum'])
    elif (isinstance(new_FVA[the_key]['maximum'],NoneType) != isinstance(reference_FVA[the_key]['maximum'],NoneType)):
        unsolved_max_one[the_key] = {}
        unsolved_max_one[the_key]['reference'] = reference_FVA[the_key]
        unsolved_max_one[the_key]['test'] = new_FVA[the_key]
    if ((not isinstance(new_FVA[the_key]['minimum'],NoneType)) & (not isinstance(reference_FVA[the_key]['minimum'],NoneType))):
        delta_dict_min[the_key] = abs(new_FVA[the_key]['minimum'] - reference_FVA[the_key]['minimum'])
    elif (isinstance(new_FVA[the_key]['minimum'],NoneType) != isinstance(reference_FVA[the_key]['minimum'],NoneType)):
        unsolved_min_one[the_key] = {}
        unsolved_min_one[the_key]['reference'] = reference_FVA[the_key]
        unsolved_min_one[the_key]['test'] = new_FVA[the_key]

# Check how many optimizations were not be solved in one case;
# This should be 0
print(len(unsolved_max_one.keys()))
print(len(unsolved_min_one.keys()))

hist(delta_dict_max.values())
hist(delta_dict_min.values())

# We can also test for required metabolites
test_metabolites = {x.id: [x] for x in gim3e_model.reactions if x.id.startswith("TMS_")}
metabolite_ko_result = gim3e.irreversible_reaction_knockout_analysis(gim3e_model,
                        the_reactions=test_metabolites,
                        solver='cplex',
                        tolerance_optimality = selected_tolerance,
                        tolerance_feasibility = selected_tolerance,
                        tolerance_barrier= 0.0001 * selected_tolerance,
                        error_reporting=None,
                        verbose = True)

# Check the results against a reference run.  
filename = "STM_LB_metabolite_ko_reference.pickle"
fp = open(data_dir + filename, "rb")
metabolite_ko_result_ref = pickle.load(fp)
fp.close()

the_cutoff = selected_tolerance
required_list = []
the_dict = metabolite_ko_result
for ko_key in the_dict.keys():
    not_required = False
    if the_dict[ko_key]['maximum_status'] not in ['failed']:
        if the_dict[ko_key]['maximum'] > the_cutoff:
            not_required = True
    if not not_required:
        required_list.append(ko_key)

# Also check against the reference dataset
required_list_ref = []
the_dict = metabolite_ko_result_ref
for ko_key in the_dict.keys():
    not_required = False
    if the_dict[ko_key]['maximum_status'] not in ['failed']:
        if the_dict[ko_key]['maximum'] > the_cutoff:
            not_required = True
    if not not_required:
        required_list_ref.append(ko_key)

# These should match
print(len(required_list))
print(len(required_list_ref))
