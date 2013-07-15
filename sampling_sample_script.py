# Make a GIM3E model from e coli core & sample
# export VERSIONER_PYTHON_PREFER_32_BIT=yes
import gim3e
import gms
import pickle
from cobra.manipulation import initialize_growth_medium
from copy import deepcopy
from types import *
from cobra.io.sbml import create_cobra_model_from_sbml_file

selected_tolerance = 1E-8
selected_growth = 0
selected_penalty = 0

# develop the sampling algorithm with E coli core as a first approach
sbml_file = 'E_coli_core_M9.xml'
cobra_model = create_cobra_model_from_sbml_file('E_coli_core_M9.xml', print_time=True)
cobra_model.reactions.get_by_id('ATPM').objective_coefficient = 1
cobra_model.optimize()
cobra_model.solution.f


gim3e.gim3e_optimize(cobra_model, objective_sense = 'maximize', 
                   the_problem = None, solver = 'cplex',  
                   error_reporting = None,
                   tolerance_optimality = selected_tolerance, 
                   tolerance_feasibility = selected_tolerance,
                   tolerance_barrier = 0.0001 * selected_tolerance,
                   tolerance_integer = 0)
cobra_model.solution.f

# Note detected metabolites that paired
# with the model are included even if
# they are blocked.  The algorithm
# automatically detects and warns
# about blocked metabolites.
# Note we can also add in virtual metabolites
# to merely monitor metabolites by formulating
# the model and not imposing flux requirements
# In a two-step process turnover flux constraints can
# then be imposed on detected metabolites.
gim3e_model, gim3e_FVA, penalty_score = gim3e.gim3e(cobra_model, expression_dict = {},
    expression_threshold = 0,
    metabolite_list = [], 
    fraction_growth = selected_growth,
    relative_penalty_bound = selected_penalty,
    solver_tolerance = selected_tolerance,
    metabolite_flux_requirement = False,
    monitor_all_cellular_metabolites = True,
    MILP_formulation = True,
    run_FVA = True, reduce_model = False,
    trim_model = False, solver = 'cplex')

sampling_object = gms.sample_container(gim3e_model)
gms.sampling_optimize(sampling_object.cobra_model_full, objective_sense = 'maximize', 
                   the_problem = None, solver = 'cplex',  
                   error_reporting = None,
                   tolerance_optimality = selected_tolerance, 
                   tolerance_feasibility = selected_tolerance,
                   tolerance_barrier = 0.0001 * selected_tolerance,
                   tolerance_integer = 0)
# make warmup points
gms.create_warmup_points(sampling_object, solver = 'cplex', solver_tolerance = 1E-8, force_vms = True, additional_forced_reactions = ['penalty'])
# get rid of redundant points
gms.reduce_warmup_points(sampling_object, solver_tolerance = selected_tolerance)
# sample
gms.achr_sampler(sampling_object, solver = "cplex", solver_tolerance = 1E-8, max_time = 60 * 60 * 24, n_points = 1000)
gms.save_sampling_object(sampling_object, "sampling_trial_core")

the_converted_results, converted_reaction_list = gms.convert_sampling_results_to_reversible(sampling_object)
the_converted_result_dict = {}
n_reactions, n_points = the_converted_results.shape
for the_row_index, the_row_id in enumerate(converted_reaction_list):
    the_converted_result_dict[the_row_id] = {}
    for the_column_index in range(0, n_points):
        the_converted_result_dict[the_row_id][the_column_index] = the_converted_results[(the_row_index, the_column_index)]
mix_frac = gms.mix_fraction(sampling_object.sampled_points, sampling_object.initial_points, fixed = sampling_object.const_ind)
print(mix_frac)

# import sys
# sys.path.append("/Users/brianschmidt/Documents/python_development")
# import dataio
# dataio.writetable(dataio.dicttotable(the_converted_result_dict, "reaction_id"),"sampling_trial_core.txt")












