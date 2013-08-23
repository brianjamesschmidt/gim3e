# Load these to help determine which solver solutions are OK
acceptable_solution_strings = ['optimal', 'MIP_optimal', 'optimal_tolerance', 'x_bound_infeasible']
# May want to include this as acceptable when running cplex
# if a small violation of the boundaries is OK: 'x_bound_infeasible'
optimal_solution_strings = ['optimal']
# Technically we can set the integer_tolerance to 0 in cplex
# This might lead to some differences in convergence from gurobi, 
# which is limited to 1e-9.
integer_tolerances = {'cplex': 0, 'gurobi': 1e-9, 'glpk': 1e-9}

def gim3e(cobra_model, expression_dict = {}, expression_threshold = 0.,
           metabolite_list = [],
           fraction_growth = 0.9,
           relative_penalty_bound = 1.,
           solver_tolerance = 1E-7,
           metabolite_flux_requirement = True,
           monitor_all_cellular_metabolites = False,
           MILP_formulation = False,
           run_FVA = True,
           reduce_model = False,
           trim_model = False,
           solver = 'cplex'):
    """ Gene Inactivity Moderated by Metabolism, Metabolomics, and Expression (GIM^3E)
    Create model from mRNA expression data (and/or potentially high coverage proteomics)
    while requiring the utilization of specific network metabolites
    
    Arguments:
     cobra_model: the model object.  note the objective should be specified in the model already
     
     expression_dict: expression data of form GENE: value
                      values are treated with linear penalty terms
                      below the threshold.  e.g. could be log2 intensities
                      or log p-values.

                      If we have p/a associated p-values,
                      can use +log p (two-tail) for the appropriate
                      directionality vs. background if absent, 0 if present,
                      with an associate threshold of log(alpha/2).  P-values
                      should be multiple testing corrected before calling
                      gim3e.

                      We can also use P/A calls, with P = 1 and A = 0
                      if we define the expression_threshold as 1
                      (not tested yet, should work)

                      If reaction ID's are given as keys, or the
                      threshold is specified as None, then these
                      (positive) penalties are applied directly to
                      the prblem (and the threshold is ignored).
                      
    expression_threshold: threshold value.  Values below threshold will be
     processed as sources of penalty.
     
    fraction_growth: fraction of the optimal objective value that must be met
     when applying constraints.  The objective should be set in the model before
     calling gim3e

    relative_penalty_bound: required fraction of the best penalty score, must
     be greater than 1.  E.g. maximum allowed penalty will be 1.05 X
     the best achievable for a given growth constraint if this is set
     to 1.05.  This is added as a model constraint.  A value of zero
     implies no limit.

    metabolite_list: list of metabolit id's to add 'turnover metabolite' constraints
     to.  To get FVA results for internal metabolites without bounds, you
     can set this to True and tm_flux_requirement to False
                 
    solver_tolerance: rough cutoff for metabolomics and solver tolerance,
      metabolites will have 1.01X tolerance 
      and optimization will have 1X tolerance.

    metabolite_flux_requirement: if True, metabolites from metabolite_list
      will be required to carry flux

    monitor_all_cellular_metabolites: if True, turnover metabolites will be added
      without flux requirements to monitor the remaining cellular metabolites

    MILP_formulation: If True, the model will be converted to a MILP.  Indicator
     variables and metabolites will be added to enforce reaction directionality.
     This eliminates flux loops for turnover metabolites but is computationally
     costly to solve.

    run_FVA: boolean, whether to run FVA for the model.  Note this must
     be set to true for reduce_model to run.
              
    reduce_model: boolean, if false we simply implement penalty as
     a constraint.  if true we further disable reactions not used
     for which data suggests are penalized and do not carry a flux.

    trim_model: boolean, whether to remove (delete) reactions that don't carry 
                 flux.                  

   Returns a tuple containing:
    new_cobra_model: the processed model.  Penalty is
      converted to a constraint, also note that the model is converted to
      irreversible to enable the constraint and monitoring of turnover
      metabolite fluxes
      
    FVA_result: the results from FVA with the new model.
     
    best_total_penalty: the best obtainable penalty measure for the specified
      constraints.


    """

    from copy import deepcopy
    from cobra.core.Reaction import Reaction
    from cobra import solvers
    import types
    continue_flag = True

    FVA_with_minimum_penalty = {}
    best_total_penalty = -1
    new_cobra_model = deepcopy(cobra_model)

    # First finalize the solver to use through the script
    solver = check_solver(solver)
    if type(solver) == types.NoneType:
        print("Cannot identify a working solver, exiting ...")
        continue_flag = False

    if continue_flag:
        # Tolerance for solver is also used as guideline for 
        # forcing metabolite fluxes
        if solver_tolerance <= 0:
            print("Warning, GIM^3E requires a positive tolerance value.")
            print("Using 1 x 10^-7")
            solver_tolerance = 1E-7
        if metabolite_flux_requirement:
            epsilon = 1.01 * solver_tolerance
        else:
            epsilon = 0
        tolerance_integer = integer_tolerances[solver]

        print("Determining the best objective value...")		
		
        gim3e_optimize(new_cobra_model, objective_sense = 'maximize', 
                   the_problem = None, solver = solver,  
                   error_reporting = None,
                   tolerance_optimality = solver_tolerance, 
                   tolerance_feasibility = solver_tolerance,
                   tolerance_barrier = 0.0001 * solver_tolerance,
                   tolerance_integer = tolerance_integer)
        if type(new_cobra_model.solution) != types.NoneType:
            best_objective = new_cobra_model.solution.f
            if new_cobra_model.solution.status not in acceptable_solution_strings:
                print("Cannot optimize the original problem, even without metabolite constraints!  Exiting...")
                continue_flag = False
            else:
                print("... OK")
        else:
            print("Cannot optimize the original problem, even without metabolite constraints!  Exiting...")
            continue_flag = False

    if continue_flag:
        if fraction_growth >= ((best_objective - solver_tolerance) / best_objective):
            fraction_growth = (best_objective - solver_tolerance) / best_objective                
            print("Fraction of growth must take numerical limitations into consideration and cannot be set too close to 1.  Resetting to " + str(fraction_growth) + ".")	
			
        # First make the model irreversible.
        # First set exchange constraints to +/- max so media exchanges are
        # are given a _reverse reaction even if it is not used.
        # This facilitates later
        # analysis on the exchanges e.g. in case we want to add
        # reactions to the media.
        print("Converting to irreversible...")
        exchange_reactions = new_cobra_model.reactions.query("EX_")
        unmodified_exchange_reactions = deepcopy(exchange_reactions)
        for the_reaction in exchange_reactions:
            the_reaction.upper_bound = 1000
            the_reaction.lower_bound = -1000
        convert_to_irreversible_with_genes(new_cobra_model, mutually_exclusive_directionality_constraint = MILP_formulation)
        for reference_reaction in unmodified_exchange_reactions:
            the_forward_reaction = new_cobra_model.reactions.get_by_id(reference_reaction.id)
            the_reverse_reaction = new_cobra_model.reactions.get_by_id(reference_reaction.id + "_reverse")        
            if reference_reaction.lower_bound < 0:
                the_reverse_reaction.lower_bound = -1 * min(reference_reaction.upper_bound, 0)
                the_reverse_reaction.upper_bound = -1 * reference_reaction.lower_bound
            else:
                the_reverse_reaction.upper_bound, the_reverse_reaction.lower_bound = (0, 0)
            if reference_reaction.upper_bound > 0:
                the_forward_reaction.upper_bound = reference_reaction.upper_bound
                the_forward_reaction.lower_bound = max(0, reference_reaction.lower_bound)
            else:
                the_forward_reaction.upper_bound, the_forward_reaction.lower_bound = (0, 0)
        print("... OK")

    if ((len(metabolite_list) > 0) & (continue_flag == True)):
        print("Adding turnover metabolites...")
        # Now add turnover metabolites for the optimization problem.
        # These are dummy metabolites that ensure flux through the
        # metabolites of interest
        add_turnover_metabolites(new_cobra_model, metabolite_list, epsilon)

        # Run FVA and verify whether any turnover metabolites are blocked.
        # First reset the TM reaction bounds to zero
        the_FVA_reactions = new_cobra_model.reactions.query('^TMS_')        
        for reaction in the_FVA_reactions:
            reaction.lower_bound = 0
        gim3e_optimize(new_cobra_model, objective_sense = 'maximize', 
                   the_problem=None, solver=solver,  
                   error_reporting=None,
                   tolerance_optimality = solver_tolerance, 
                   tolerance_feasibility = solver_tolerance,
                   tolerance_barrier=0.0001 * solver_tolerance,
                   tolerance_integer = tolerance_integer)

        # Now that we have a MILP, use the gim3e_optimize function,
        # which controls additional solver parameters.
        gim3e_optimize(new_cobra_model,objective_sense = 'maximize', 
            the_problem=None, solver=solver,  
            error_reporting=None,
            tolerance_optimality = solver_tolerance, 
            tolerance_feasibility = solver_tolerance,
            tolerance_barrier=0.0001 * solver_tolerance,
            tolerance_integer = tolerance_integer)
        best_objective_with_mets = (new_cobra_model.solution.f)
		
        # This flag should never show up, but here just in case issues arise that require troubleshooting
        if new_cobra_model.solution.status not in acceptable_solution_strings:
            print("Cannot optimize the problem with turnover metabolites, even with no minimum flux requirement!  Exiting...")
            continue_flag = False
        else:
            # We only need to check all of the turnover metabolites if
            # we are constraining with "nonzero" flux
            if epsilon > 0:
                test_turnover_metabolite_constraints(new_cobra_model,
                    epsilon = epsilon,
                    fraction_growth = fraction_growth,
                    objective_sense = 'maximize',
                    the_FVA_reactions = the_FVA_reactions,
                    solver = solver,
                    tolerance_optimality = solver_tolerance,
                    tolerance_feasibility = solver_tolerance,
                    tolerance_barrier = 0.0001 * solver_tolerance,
                    tolerance_integer = tolerance_integer,
                    error_reporting = None,
                    number_of_processes = 1,
                    FVA_verbose = False)
        print("... OK, done creating and verifying turnover metabolites (individually).")
    elif continue_flag == True:
        best_objective_with_mets = best_objective
		
    if (monitor_all_cellular_metabolites & (continue_flag == True)):
        cellular_metabolite_list = []
        for the_metabolite in new_cobra_model.metabolites:
            if ((the_metabolite.compartment == 'c') | (the_metabolite.compartment == 'p')):
                if the_metabolite.id not in metabolite_list:
                    cellular_metabolite_list.append(the_metabolite.id)
        add_turnover_metabolites(new_cobra_model, cellular_metabolite_list, 0.)

    if (continue_flag == True):
        # Now evaluate the penalties
        penalties = evaluate_penalties(cobra_model, new_cobra_model, expression_dict, expression_threshold)
        # Change the existing objective to a constraint.
        if (best_objective_with_mets < (fraction_growth * best_objective)):
            constraint_test_model = convert_objective_to_constraint(new_cobra_model,
                objective_sense = 'maximize', 
                fraction_of_optimum = fraction_growth * best_objective,
                copy_model = True,
                bound_best_optimum = False,
                new_reaction_name = "test_objective",
                solver=solver,
                tolerance_optimality = solver_tolerance,
                tolerance_feasibility = solver_tolerance,
                tolerance_barrier = 0.0001 * solver_tolerance,
                tolerance_integer = tolerance_integer)
            constraint_test_model.reactions.get_by_id("test_objective").objective_coefficient = 1.
            if (constraint_test_model.solution.status != 'optimal'):
                fraction_growth = best_objective_with_mets / best_objective
                print("Warning: cannot operate at desired fraction of optimum.  Reducing to: "+ str(fraction_growth) +" of value with no metabolite constraints.")
            elif (constraint_test_model.solution.f < fraction_growth * best_objective):
                # Check about relaxing constraint; for practical reasons it does not make
                # sense to relax it less than the solver tolerance
                # Also, GIM^3E does not perform well with fraction_growth == 1
                fraction_growth = min(best_objective_with_mets / best_objective, (best_objective - solver_tolerance) / best_objective)
                print("Warning: cannot operate at desired fraction of optimum.  Reducing to: "+ str(fraction_growth) +" of value with no metabolite constraints.")


        # Use the reaction penalties as an objective.
        old_objective = {}      
        for reaction in new_cobra_model.reactions:
            if reaction.objective_coefficient != 0:
                old_objective[reaction.id] = reaction.objective_coefficient           
        if (len(old_objective.keys()) < 2):
            the_growth_objective_id = old_objective.keys()[0]
            old_objective_reaction = new_cobra_model.reactions.get_by_id(the_growth_objective_id)
            the_growth_objective_coefficient = old_objective_reaction.objective_coefficient
            old_objective_reaction.objective_coefficient = 0
            old_objective_reaction.lower_bound = fraction_growth * best_objective
        else:
            new_cobra_model = convert_objective_to_constraint(new_cobra_model, objective_sense = 'maximize', 
                fraction_of_optimum = fraction_growth, copy_model = True,
                new_reaction_name = "objective",
                solver=solver,
                tolerance_optimality = solver_tolerance,
                tolerance_feasibility = solver_tolerance,
                tolerance_barrier = 0.0001 * solver_tolerance,
                tolerance_integer = tolerance_integer,
                bound_best_optimum = False)
            the_growth_objective_id = new_cobra_model.reactions[(len(new_cobra_model.reactions))-1].id
            the_growth_objective_coefficient = 1.
        for reaction in new_cobra_model.reactions:
            if penalties[reaction.id] > 0:
                reaction.objective_coefficient = penalties[reaction.id]       
        # Also optimize the model for the the penalty score.
        # Penalties are positive so this becomes a
        # minimization problem.

        if (sum(penalties.values()) > 0):
            # It may not be reasonable to use very "tight"
            # tolerances when optimizing e.g. the
            # penalty, it can grow to be very large
            gim3e_optimize(new_cobra_model, objective_sense = 'minimize', 
                   the_problem = None, solver = solver,  
                   error_reporting = None,
                   tolerance_optimality = solver_tolerance, 
                   tolerance_feasibility = solver_tolerance,
                   tolerance_barrier = 0.0001 * solver_tolerance,
                   tolerance_integer = tolerance_integer)
            if type(new_cobra_model.solution) != types.NoneType:
                best_total_penalty = new_cobra_model.solution.f
                if new_cobra_model.solution.status not in acceptable_solution_strings:	
                    print("Failed to find a solution when minimizing the penalties, exiting...")
                    continue_flag = False
            else:
                print("Failed to find a solution when minimizing the penalties, exiting...")
                continue_flag = False
        else:
            best_total_penalty = 0
            print("No microarray data given, not imposing constraint based on intensity data.")
            # In the case of no penalties we should also prematurely reset the
            # objective to the original and overwrite the relative_penalty_bound
            # as zero.  Thus when we run FVA, we will have identified a
            # solution for hot-starting and won't have errors in FVA from the
            # penalty bounds since growth reaction bounds are added a constraint.
            new_cobra_model.reactions.get_by_id(the_growth_objective_id).objective_coefficient = 1
            relative_penalty_bound = 0
    else:
        best_total_penalty = 0

    if (continue_flag == True):
        if relative_penalty_bound <= 0:
            test_bound = 0
        else:
            # Invert to be compatible with the irreversible_flux_variability_analysis
            test_bound = 1 / relative_penalty_bound
        if run_FVA:
            print("About to run FVA.")
            # If we have made an irreversible model, we don't want to run FVA on the
            # indicator reactions
            the_reactions = [x.id for x in new_cobra_model.reactions if (not x.id.startswith("IRRMILP_"))]
            # Track progress, this may take a while
            FVA_verbose = True
            # Run FVA and do not allow the penalty score to degrade below test_bound
            FVA_with_minimum_penalty = irreversible_flux_variability_analysis(new_cobra_model,
                fraction_of_optimum = test_bound,
                objective_sense='minimize',
                the_reactions = the_reactions,
                solver=solver,
                tolerance_optimality = solver_tolerance,
                tolerance_feasibility = solver_tolerance,
                tolerance_barrier = 0.0001*solver_tolerance,
                tolerance_integer = tolerance_integer,
                error_reporting = None,
                number_of_processes = 1,
                verbose = FVA_verbose)
        else:
            FVA_with_minimum_penalty = {}

        # We also want to convert the penalty to a constraint and the 
        # growth/ original back to an objective before returning the model
        if (len(penalties) > 0) & (sum(penalties.values()) > 0):
            new_cobra_model = convert_objective_to_constraint(new_cobra_model, objective_sense = 'minimize', 
                fraction_of_optimum = test_bound,
                copy_model = False,
                new_reaction_name = "penalty",
                solver = solver,
                tolerance_optimality = solver_tolerance,
                tolerance_feasibility = solver_tolerance,
                tolerance_barrier = 0.0001 * solver_tolerance,
                tolerance_integer = tolerance_integer,
                bound_best_optimum = False)        
        new_cobra_model.reactions.get_by_id(the_growth_objective_id).objective_coefficient = the_growth_objective_coefficient
    else:
        # No FVA results
        FVA_with_minimum_penalty = {}

    if ((continue_flag == True) & (run_FVA == True) & (reduce_model == True)):
        # We set the boundaries on reactions not included in the solution to zero
        # if expression data suggests they are absent.  Couple the forward
        # and reverse reactions to the same fate.
        # Here, we only turn off reactions if there is no growth or
        # penalty, even if we allowed a bit more slack in the FVA solutions.
        # Not pressure tested yet.
        inactivation_candidates = [x for x in penalties.keys() if penalties[x] > 0]
        
        # We don't want to turn off a reaction if it is required the current optimal growth
        # In theory this step could result in a worsening of the lowest penalty
        # within the  previous FVA bounds at the optimal as it is currently formulated
        # Use the solver tolerance as a guidline for evaluation, if reactions are carrying flux
        # greater than the solver tolerance we don't need to evaluate for inactivation.
        constrained_reaction_list = []
        
        new_cobra_model, inactivated_reactions = turn_off_unused_reactions(new_cobra_model,
            FVA_dict = FVA_with_minimum_penalty,
            epsilon = solver_tolerance,
            solver = solver,
            tolerance_optimality = solver_tolerance,
            tolerance_feasibility = solver_tolerance,
            tolerance_barrier = 0.0001*solver_tolerance,
            tolerance_integer = tolerance_integer,
            candidate_reaction_id_list = inactivation_candidates,
            # Don't allow the optimum to decrease due to
            # the removal of reactions
            fraction_of_optimum = 1.,
            # Make sure the penalty is not impacted
            constrained_reaction_list = constrained_reaction_list,
            verbose = False)
        if trim_model:
            new_cobra_model, FVA_with_minimum_penalty = remove_model_reactions(new_cobra_model,
                inactivated_reactions, FVA_with_minimum_penalty)

    # Finally, it is nicer to return a model object with the solution
    if continue_flag:
        gim3e_optimize(new_cobra_model, objective_sense = 'maximize', 
                   the_problem = None, solver = solver,  
                   error_reporting = None,
                   tolerance_optimality = solver_tolerance, 
                   tolerance_feasibility = solver_tolerance,
                   tolerance_barrier = 0.0001 * solver_tolerance,
                   tolerance_integer = tolerance_integer)

    return (new_cobra_model, FVA_with_minimum_penalty, best_total_penalty)


def add_turnover_metabolites(cobra_model, metabolite_id_list, epsilon):
    """ NOTE: Model must first be converted to irreversible!
    This entry creates a corresponding turnover metabolite
    that ensures flux through the metabolite of interest.

    Arguments:
     cobra_model: the model to be updated.
     metabolite_id_list: list of model metabolites for
                         which to add a turnover metabolite.
     epsilon: minimal flux to force through turnover metabolites.
      recommend 1.01 X solver_tolerance    

    
    """
    from cobra.core.Reaction import Reaction
    from cobra.core.Metabolite import Metabolite
    from numpy import abs

    # Set the minimum flux for metabolites equal to some factor larger than the solver's tolerance
    the_min_flux = epsilon

    turnover_metabolites = []
    sink_reactions = []
    for metabolite_id in metabolite_id_list:
        v_metabolite = Metabolite("TM_" + metabolite_id)
        # Now for reactions.  We include all reactions 
        # that create or consume the real metabolite.
        # These reactions therefore also drive creation of the
        # turnover metabolite, and we need to add a reaction
        # that is the sink.  By constraining this single
        # sink reaction, we ensure flux through the real reactions
        # to create the real metabolite.
        r_metabolite = cobra_model.metabolites.get_by_id(metabolite_id)
        sum_abs_source_reaction_bounds = 0

        the_reaction_id_list = [x.id for x in r_metabolite.get_reaction()]
        for the_reaction_id in the_reaction_id_list:
            the_reaction = cobra_model.reactions.get_by_id(the_reaction_id)
            coefficient = abs(the_reaction.get_coefficient(r_metabolite))
            the_reaction._metabolites[v_metabolite] = coefficient
            v_metabolite._reaction.add(the_reaction)
            # The model should be irreversible, and the upper bound
            # should be greater than the lower bound and positive.
            # Use 1000 for each reaction to be conservative.
            # E.g. might flip reactions back on later, don't
            # use the current upper bound for the reactions.
            sum_abs_source_reaction_bounds += 1000.     

        # Add the sink reaction for the turnover metabolite
        # Since both creation and consumption of
        # the real metabolite drive the turnover,
        # we require 2 units of metabolite per
        # 1 unit of flux sink so this matches
        # the turnover through the real metabolite.
        sink_reaction = Reaction("TMS_" + metabolite_id)
        sink_reaction.add_metabolites({v_metabolite:-2})

        # Ensure a positive flux through the sink.
        # and maximum equal to maximal needed to
        # sustain reactions at their maximum.
        sink_reaction.lower_bound = the_min_flux
        sink_reaction.upper_bound = sum_abs_source_reaction_bounds / 2.
        v_metabolite._reaction.add(sink_reaction)        
        turnover_metabolites.append(v_metabolite)
        sink_reactions.append(sink_reaction)
    cobra_model.add_metabolites(turnover_metabolites)
    cobra_model.add_reactions(sink_reactions)


def convert_to_irreversible_with_genes(cobra_model, mutually_exclusive_directionality_constraint = False):
    """Will break all of the reversible reactions into two separate irreversible
     reactions with different directions.  This function call modified from
     a version in the core cobra to facilitate the MILP formulation and
     include gene_reaction_rules with the reverse reaction

     Arguments:
      cobra_model: A model object which will be modified in place.
      mutually_exclusive_directionality_constraint: Boolean.  If True, turnover 
       reactions are constructed to serve as MILP constraints to prevent loops.
      
     Returns:
      None, cobra_model is modified in place
    
    NOTE: This function has been modified from the manipulate module

    
    """
    reactions_to_add = []
    from cobra.core.Reaction import Reaction
    from cobra.core import Metabolite    
    for reaction in cobra_model.reactions:
        # Potential artifact because a reaction might run backwards naturally
        # and this would result in adding an empty reaction to the
        # model in addition to the reverse reaction.
        if reaction.lower_bound < 0:
            #reverse_reaction = Reaction(reaction.id + "_reverse")
            reverse_reaction = reaction.copy()
            reverse_reaction.id = reaction.id + "_reverse"
            reverse_reaction.lower_bound = 0
            reverse_reaction.upper_bound = reaction.lower_bound * -1.
            reaction.lower_bound = 0
            # Make the directions aware of each other
            reaction.reflection = reverse_reaction
            reverse_reaction.reflection = reaction
            reaction.reversibility = 0
            reverse_reaction.reversibility = 0
            reaction_dict = {}
            current_metabolites = [x for x in (reaction.get_products() + reaction.get_reactants())]
            for the_metabolite in current_metabolites:
                reaction_dict[the_metabolite] = -2 * reaction.get_coefficient(the_metabolite.id)
            reverse_reaction.add_metabolites(reaction_dict)
            reactions_to_add.append(reverse_reaction)
            # Also: GPRs should already copy
            # reverse_reaction.gene_reaction_rule = reaction.gene_reaction_rule
            # reverse_reaction._genes = reaction._genes
            if mutually_exclusive_directionality_constraint:
                # A continuous reaction bounded by 0., 1.
                # Serves as a source for the indicator metabolites
                tmp_source = Reaction('IRRMILP_direction_constraint_source_for_%s_and_%s'
                                                           %(reaction.id,
                                                             reverse_reaction.id))
                tmp_source.upper_bound = 1.
                tmp_source.lower_bound = 0.
                # The reverse indicator reaction is
                # an integer-valued reaction bounded by 0,1
                # that activates flux to the reverse reaction
                # and deactivates the forward reaction only when it is
                # turned on to 1
                tmp_indicator = Reaction('IRRMILP_reverse_indicator_for_%s_and_%s'
                                                           %(reaction.id,
                                                             reverse_reaction.id))
                tmp_indicator.upper_bound = 1
                tmp_indicator.lower_bound = 0
                tmp_indicator.variable_kind = 'integer'                    
                flux_constraint_forward = Metabolite(id = 
                     'IRRMILP_direction_constraint_for_%s'%reaction.id)
                flux_constraint_reverse = Metabolite(id = 
                     'IRRMILP_direction_constraint_for_%s'%reverse_reaction.id)
                flux_constraint_reverse._constraint_sense = 'G'
                flux_constraint_reverse._bound = 0.

                tmp_source.add_metabolites({flux_constraint_forward: 1})

                tmp_indicator.add_metabolites({flux_constraint_forward: -1,
                                      flux_constraint_reverse: 1})
                if reaction.upper_bound != 0:
                        reaction.add_metabolites({flux_constraint_forward: -1./reaction.upper_bound})
                else:
                    # could put 1.01 X the tolerance here,
                    # This is arbitrary.  Use 0.001
                    # since 1000 is a typical upper bound
                    reaction.add_metabolites({flux_constraint_forward: -0.001})
                if reverse_reaction.upper_bound != 0:
                    reverse_reaction.add_metabolites({flux_constraint_reverse: -1./reverse_reaction.upper_bound})
                else:
                    reverse_reaction.add_metabolites({flux_constraint_reverse: -0.001})
                reactions_to_add.append(tmp_indicator)
                reactions_to_add.append(tmp_source)
    cobra_model.add_reactions(reactions_to_add)


def remove_model_reactions(cobra_model, the_reactions_to_remove = [], FVA_result_dict = {}):
    """ Trim reactions from a model from a model.
    After calling turn_off_unused_reactions, this function can be used to
    delete the selected reactions and update the FVA_result_dict.
    Note if the reactions are not zero-flux incongruency may result
    between the model and result dict

    Arguments:
     cobra_model: model object, modified in place
     FVA_result_dict: Optional FVA result, keys are reaction id's.
      not used for decision making, 
     the_reactions_to_remove: a list of reaction id's to remove

    Output:
     cobra_model: modified in place, too.
     FVA_result_dict: modified in place, too,
      pop entries, won't re-run.


    """
    from copy import deepcopy
    from cobra.core.Reaction import Reaction
    
    approved_reactions = []
    for reaction_id in the_reactions_to_remove:
        # Remove zero-flux reactions.  Just
        # Do a check to make sure we just add one reaction
        # from each pair.
        the_reaction = cobra_model.reactions.get_by_id(reaction_id)
        if the_reaction not in approved_reactions:
            if 'reflection' in dir(the_reaction):
                if type(the_reaction.reflection) == Reaction:
                    if the_reaction.reflection.id not in the_reactions_to_remove:
                        approved_reactions.append(the_reaction)
                        # Just add one from each pair
                else:
                    approved_reactions.append(the_reaction)
            else:
                approved_reactions.append(the_reaction)
    # Now trim the model
    IRRMILP_id_list = [x.id for x in cobra_model.reactions if x.id.startswith("IRRMILP_")]
    for the_reaction in approved_reactions:
        if the_reaction.id in FVA_result_dict.keys():
            FVA_result_dict.pop(the_reaction.id)
        if 'reflection' in dir(the_reaction):
            if type(the_reaction.reflection) == Reaction:
                if the_reaction.reflection.id in FVA_result_dict.keys():
                    FVA_result_dict.pop(the_reaction.reflection.id)
                comb_1 = 'IRRMILP_direction_constraint_source_for_'+the_reaction.id+'_and_'+the_reaction.reflection.id
                comb_2 = 'IRRMILP_direction_constraint_source_for_'+the_reaction.reflection.id+'_and_'+the_reaction.id
                if (comb_1 in IRRMILP_id_list):
                    IRRMILP_direction_constraint_source_reaction = cobra_model.reactions.get_by_id(comb_1)
                    IRRMILP_direction_constraint_source_reaction.remove_from_model(cobra_model)
                elif (comb_2 in IRRMILP_id_list):
                    IRRMILP_direction_constraint_source_reaction = cobra_model.reactions.get_by_id(comb_2)
                    IRRMILP_direction_constraint_source_reaction.remove_from_model(cobra_model)
                comb_1 = 'IRRMILP_reverse_indicator_for_'+the_reaction.id+'_and_'+the_reaction.reflection.id
                comb_2 = 'IRRMILP_reverse_indicator_for_'+the_reaction.reflection.id+'_and_'+the_reaction.id
                if (comb_1 in IRRMILP_id_list):
                    IRRMILP_reverse_indicator_for_reaction = cobra_model.reactions.get_by_id(comb_1)
                    IRRMILP_reverse_indicator_for_reaction.remove_from_model(cobra_model)
                elif (comb_2 in IRRMILP_id_list):
                    IRRMILP_reverse_indicator_for_reaction = cobra_model.reactions.get_by_id(comb_2)
                    IRRMILP_reverse_indicator_for_reaction.remove_from_model(cobra_model)
                the_reaction.reflection.remove_from_model(cobra_model)

            the_reaction.remove_from_model(cobra_model)
            # print(the_reaction.id)

      
    # Remove orphaned metabolites
    for metabolite in cobra_model.metabolites:
        if len(metabolite.get_reaction()) == 0:
            metabolite.remove_from_model(cobra_model)
    # Remove orphan genes
    for gene in cobra_model.genes:
        if len(gene.get_reaction()) == 0:
            gene.remove_from_model(cobra_model)

    return (cobra_model, FVA_result_dict)


def irreversible_flux_variability_analysis(cobra_model, **kwargs):
    """
    Perform flux variability analysis on a model that has been converted to irreversible.
    
    Arguments:
     cobra_model: Preferably, an objecive reaction should be identified.  If not, we will
      just pick one to get a solution for hot-starting.  

    kwargs:
     fraction_of_optimum: if an objective is present, it will be converted to a constraint and
      the fraction_of_optimum used to set the boundary.

     objective_sense: 'minimize' or 'maximize'

     the_reactions: a list of reactions to perform FVA with.  If not declared,
      all model reactions will be tested.  If a dict is supplied with
      reaction ids as keys, FVA will be run using the supplied values
      as upper bounds.  The upper limit on the flux for the
      the target reactions will be changed to this value.

     Next arguments are solver parameters, described previously:
      solver, tolerance_optimality, tolerance_feasibility, tolerance_barrier, tolerance_integer, error_reporting, number_of_processes

     cold_start_list: reaction id's to force a cold start.  A safety feature,
      very rarely the solver has frozen when hot starting.

     verbose: if true, output the current reaction.  Useful for long
      e.g. MILP computations as well as troubleshooting reactions.

    Returns:
     variability_dict: a dictionary with the reaction id's as top-level keys,
      FVA ranges, and solver status

    
    """
    from numpy import hstack, zeros, array
    from scipy import sparse
    from copy import deepcopy
    from cobra.core.Metabolite import Metabolite
    from cobra.core.Reaction import Reaction
    from time import time
    import types
    from math import floor

    # Set up
    if 'fraction_of_optimum' in kwargs:
        fraction_of_optimum = kwargs['fraction_of_optimum'] 
    else:
        fraction_of_optimum = 1.

    if 'objective_sense' in kwargs:
        objective_sense = kwargs['objective_sense'] 
    else:
        objective_sense = 'maximize'

    if 'the_reactions' in kwargs:
        the_reactions = kwargs['the_reactions'] 
    else:
        the_reactions = None

    if 'error_reporting' in kwargs:
        error_reporting = kwargs['error_reporting'] 
    else:
        error_reporting = None

    if 'number_of_processes' in kwargs:
        number_of_processes = kwargs['number_of_processes'] 
    else:
        number_of_processes = 1

    solver, tolerance_optimality, tolerance_feasibility, tolerance_barrier, tolerance_integer = get_solver_parameters(**kwargs)

    if 'cold_start_list' in kwargs:
        cold_start_list = kwargs['cold_start_list'] 
    else:
        cold_start_list = []

    if 'verbose' in kwargs:
        verbose = kwargs['verbose'] 
    else:
        verbose = False

    the_bounds = None
    continue_flag = True

    # INITIALIZATION
    # Need to copy the model because we're updating
    # reactions.
    cobra_model = cobra_model.copy()
    if not the_reactions:
        the_reactions = cobra_model.reactions

    elif type(the_reactions) == dict:
        # Break this into two lists, one of reactions and one of bounds.
        the_bounds = []
        the_reactions_temp = []
        for the_reaction in the_reactions:
            if hasattr(the_reaction, 'id'):
                the_reactions_temp.append(the_reaction.id)
            else:
                the_reactions_temp.append(the_reaction)
            the_bounds.append(the_reactions[the_reaction])
        the_reactions = map(cobra_model.reactions.get_by_id, the_reactions_temp)

    else:
        if hasattr(the_reactions[0], 'id'):
            the_reactions = [x.id for x in the_reactions]
        the_reactions = map(cobra_model.reactions.get_by_id, the_reactions)

    # We will need to keep these in memory.
    if solver == 'cplex':
        from cobra.solvers.cplex_solver import get_status
    elif solver == 'gurobi':
        from cobra.solvers.gurobi_solver import get_status
    elif solver == 'glpk':
        from cobra.solvers.glpk_solver import get_status
    
    variability_dict = {}
    # Not yet supported
    if number_of_processes > 1:
        print("This FVA script developed for irreversible models is not set up for multiple processors.  Exiting...")
        continue_flag = False

    # Reformulate the optimization problem.
    else:
        # If an objective is detected, convert it to
        # a constraint for FVA to use the fr of optimum
        nonzero_objective_list = [x.id for x in cobra_model.reactions if abs(x.objective_coefficient) > 0]
        nonzero_objective_coefficient = 1
        if (len(nonzero_objective_list) > 0):
            if (len(nonzero_objective_list) > 1):            
                # Numerically it is safer to leave the optimum
                # implicitly constrained rather than adding a constraint
                # that essentially repeats this.  This was observed to make
                # a difference for rare MILP optimizations that
                # could not be solved in cplex when the redundant constraint
                # was added in.    
                convert_objective_to_constraint(cobra_model,
                    objective_sense = objective_sense, 
                    fraction_of_optimum = fraction_of_optimum,
                    copy_model = False,
                    new_reaction_name = "objective",
                    solver=solver,
                    tolerance_optimality = tolerance_optimality,
                    tolerance_feasibility = tolerance_feasibility,
                    tolerance_barrier = tolerance_barrier,
                    tolerance_integer = tolerance_integer,
                    bound_best_optimum = False)
                objective_reaction = cobra_model.reactions.get_by_id('objective')
                objective_reaction.objective_coefficient = nonzero_objective_coefficient
                nonzero_objective_list = [objective_reaction.id]
            else:
                the_output = add_sense_aware_bounds(cobra_model,
                    objective_sense = objective_sense, 
                    fraction_of_optimum = fraction_of_optimum,
                    solver=solver,
                    tolerance_optimality = tolerance_optimality,
                    tolerance_feasibility = tolerance_feasibility,
                    tolerance_barrier = tolerance_barrier,
                    tolerance_integer = tolerance_integer,
                    bound_best_optimum = False)
                if the_output != None:
					objective_reaction = cobra_model.reactions.get_by_id(nonzero_objective_list[0])
					nonzero_objective_coefficient = objective_reaction.objective_coefficient
                else:
					print("Unable to add a minimum bound for the objective reaction in FVA. Exiting ...")
					continue_flag = False
            best_solution = gim3e_optimize(cobra_model,
                solver=solver,
                objective_sense=objective_sense,
                tolerance_optimality=tolerance_optimality,
                tolerance_feasibility=tolerance_feasibility,
                tolerance_barrier = tolerance_barrier,
                tolerance_integer = tolerance_integer)            
        else:        
            # if so there's no objective
            # set to get a hot-start
            # solution for FVA runs           
            candidates = cobra_model.reactions.query("biomass")
            objective_reaction = None
            for test_reaction in candidates:
                if ((test_reaction.upper_bound != 0) | (test_reaction.lower_bound != 0)):
                    objective_reaction = test_reaction
            if objective_reaction == None:
                print('Unable to identify an objective for hot-starting FVA solutions.  Exiting.')
                continue_flag = False
                best_solution = None
            else:
                objective_reaction.objective_coefficient = 1.
                nonzero_objective_coefficient = 1
                nonzero_objective_list = [objective_reaction.id]
                best_solution = gim3e_optimize(cobra_model,
                    solver = solver,
                    objective_sense = objective_sense,
                    new_objective = objective_reaction,
                    tolerance_optimality = tolerance_optimality,
                    tolerance_feasibility = tolerance_feasibility,
                    tolerance_barrier = tolerance_barrier,
                    tolerance_integer = tolerance_integer)
        
        if continue_flag:
            # Back up the optimal solution to access later
            optimized_model = deepcopy(cobra_model)
            objective_reaction.objective_coefficient = 0.
            if verbose:
                start_time = time()
            # FVA: STEP THROUGH REACTIONS
            the_sense_dict = {'maximize': 'maximum',
                          'minimize': 'minimum'}
            for the_reaction in the_reactions:
                # This flag makes sure we try to find a solution with the best basis
                move_to_next = False
                # We try to start with the best solution
                # but will have to start froms scratch if we violate it
                # by altering the reverse reaction bounds.
                the_solution = best_solution                
                the_reaction = cobra_model.reactions.get_by_id(the_reaction.id)
                tmp_dict = {}

                if the_bounds != None:
                    reaction_upper_bound = the_reaction.upper_bound
                    the_reaction.upper_bound = the_bounds[the_reactions.index(the_reaction)]
                # e.g. To explore variability in the forward reaction we constrain
                # the reverse reaction to zero.
                if 'reflection' in dir(the_reaction):
                    if type(the_reaction.reflection) == Reaction:
                        reflection_id = the_reaction.reflection.id
                        reflection = cobra_model.reactions.get_by_id(reflection_id)
                        # Note there is a potential artifact: if the reverse
                        # reaction must have nonzero flux. As a patch,
                        # the solution will be reset to infeasible later if we violate this
                        reflection_lower_bound = cobra_model.reactions.get_by_id(reflection_id).lower_bound
                        cobra_model.reactions.get_by_id(reflection_id).lower_bound = 0
                        reflection_upper_bound = cobra_model.reactions.get_by_id(reflection_id).upper_bound
                        cobra_model.reactions.get_by_id(reflection_id).upper_bound = 0
                        if (optimized_model.solution.x_dict[reflection_id] != 0.):
                            # If the old solution now violates a
                            # constraint we will need to start
                            # the solver from scratch.  E.g. have
                            # seen an artifact such that GLPK will pass
                            # off a solution in violation with an acceptable
                            # flag if we start in violation.
                            the_solution = None
                if the_reaction.id in cold_start_list:
                    the_solution = None
                # Old objectives should have been zero'd out and
                # incorporated into production of the objective metabolite
                the_reaction.objective_coefficient = 1.
                
                for the_sense, the_description in the_sense_dict.iteritems():
                    move_to_next = False
                    while not move_to_next:
                        gim3e_optimize(cobra_model,
                            solver = solver,
                            objective_sense = the_sense,
                            # Need to include new_objective
                            # here to update the_problem
                            # when hot-starting.
                            new_objective = the_reaction,
                            tolerance_optimality = tolerance_optimality,
                            tolerance_feasibility = tolerance_feasibility,
                            tolerance_barrier = tolerance_barrier,
                            tolerance_integer = tolerance_integer,
                            the_problem = the_solution,
                            error_reporting = error_reporting)                                
                        if type(cobra_model.solution) != types.NoneType:
                            tmp_dict[the_description] = cobra_model.solution.f
                            tmp_dict[str(the_description)+"_"+"status"] = cobra_model.solution.status
                            move_to_next = True
                        else:
                            # Give this one more try if we haven't reset the basis
                            if the_solution != None:
                                the_solution = None
                                move_to_next = False
                            else:
                                move_to_next = True
                                tmp_dict[the_description] = None
                                tmp_dict[str(the_description)+"_"+"status"] = 'infeasible'

                # Additional check for solutions
                # If one max/min optimize was successful but not the other,
                # we can retry the failed one using the solution from the
                # successful one as a basis
                n_bad_attempts = 0
                for the_sense, the_description in the_sense_dict.iteritems():
                    if tmp_dict[str(the_description)+"_"+"status"] not in acceptable_solution_strings:
                        n_bad_attempts += 1
                if n_bad_attempts == 1:
                    for the_sense, the_description in the_sense_dict.iteritems():
                        if tmp_dict[str(the_description)+"_"+"status"] not in acceptable_solution_strings:
                            failed_sense = the_sense
                            failed_description = the_description
                        else:
                            successful_sense = the_sense
                    successful_solution = gim3e_optimize(cobra_model,
                            solver = solver,
                            objective_sense = successful_sense,
                            # Need to include new_objective
                            # here to update the_problem
                            # when hot-starting.
                            new_objective = the_reaction,
                            tolerance_optimality = tolerance_optimality,
                            tolerance_feasibility = tolerance_feasibility,
                            tolerance_barrier = tolerance_barrier,
                            tolerance_integer = tolerance_integer,
                            the_problem = the_solution,
                            error_reporting = error_reporting)
                    if type(cobra_model.solution) != types.NoneType:
                        if cobra_model.solution.status in acceptable_solution_strings:
                            gim3e_optimize(cobra_model,
                                                solver = solver,
                                                objective_sense = failed_sense,
                                                # Need to include new_objective
                                                # here to update the_problem
                                                # when hot-starting.
                                                new_objective = the_reaction,
                                                tolerance_optimality = tolerance_optimality,
                                                tolerance_feasibility = tolerance_feasibility,
                                                tolerance_barrier = tolerance_barrier,
                                                tolerance_integer = tolerance_integer,
                                                the_problem = successful_solution,
                                    error_reporting = error_reporting)
                            if type(cobra_model.solution) != types.NoneType:
                                tmp_dict[failed_description] = cobra_model.solution.f
                                tmp_dict[str(failed_description)+"_"+"status"] = cobra_model.solution.status

                # Need to perform a check to deal with 
                # numerical issues due to sub-tolerance artifacts.
                tmp_dict = check_bounds_consistency(tmp_dict, the_reaction)
                the_reaction.objective_coefficient = 0
                if the_bounds != None:
                    the_reaction.upper_bound = reaction_upper_bound
                    # The lower bound isn't altered in the_bounds
                # Now reset the reaction bounds for the reverse reaction.
                if 'reflection' in dir(the_reaction):
                    if type(the_reaction.reflection) == Reaction:
                        cobra_model.reactions.get_by_id(reflection_id).lower_bound = reflection_lower_bound
                        cobra_model.reactions.get_by_id(reflection_id).upper_bound = reflection_upper_bound
                        # Set the solution as infeasible of we have violated a requirement for
                        # the reverse reaction to carry flux; don't expect to
                        # run into this situation frequently
                        if ((reflection_lower_bound > 0) & (reflection_upper_bound > 0)):
                            for the_sense, the_description in the_sense_dict.iteritems():
                                tmp_dict[the_description] = None
                                # Solution is infeasible if it violates a boundary
                                tmp_dict[the_description + "_" + "status"] = 'infeasible'                                
                            
                variability_dict[the_reaction.id] = tmp_dict
                # Check the status for the best_solution to make sure
                # it is still valid.  If not, need to reset it.
                best_solution_status = get_status(best_solution)
                if (best_solution_status not in acceptable_solution_strings):
                    # Reset best_solution for future tries.  Note this list was converted to 1 element previously.
                    the_objective_reaction = nonzero_objective_list[0]
                    cobra_model.reactions.get_by_id(the_objective_reaction).objective_coefficient = nonzero_objective_coefficient
                    best_solution = gim3e_optimize(cobra_model,
                                                            solver = solver,
                                                            objective_sense = objective_sense,
                                                            new_objective = nonzero_objective_list[0],
                                                            tolerance_optimality = tolerance_optimality,
                                                            tolerance_feasibility = tolerance_feasibility,
                                                            tolerance_barrier = tolerance_barrier,
                                tolerance_integer = tolerance_integer)
                    cobra_model.reactions.get_by_id(the_objective_reaction).objective_coefficient = 0
                            
                if verbose:
                    next_reaction_id = "...none.  All done"
                    the_index = the_reactions.index(the_reaction.id)
                    pass_s = time() - start_time
                    remaining_s = pass_s * (len(the_reactions) - (the_index + 1)) / (the_index + 1)
                    remaining_h = floor((remaining_s)/3600)
                    remaining_m = floor(((remaining_s)/3600 - remaining_h) * 60)                    
                    if (the_index + 1) < len(the_reactions):
                        next_reaction_id = the_reactions[the_index + 1].id
                    print("Completed "+ str(the_index + 1 ) + " of " + str(len(the_reactions)) + ", about to try " + next_reaction_id + ".  El: %0.0f s.  R: %0.0f hr %0.0f min." % (pass_s, remaining_h, remaining_m))                        

    return variability_dict


def convert_FVA_dict_to_reversible(irreversible_FVA_dict, cobra_model, enforce_consistency = True):
    """ Take a dict returned by irreversible_FVA and convert it to the corresponding
    reversible result where reversible reactions have been consolidated back into one.
    
    Arguments:
     irreversible_FVA_dict: results of irreversible FVA
     cobra_model: the corresponding cobra_model that contains
                  data on the reversible reaction in
     enforce_consistency: whether to enforce numerical consistency
      with respect to ordering of the min and max and the bounds.
      It may be desirable to set to false in cases where the bounds have
      been manually over-ridden, eg with FVA.


    Returns:
     consolidated_dict: the dict with the reversible reactions consolidated
     
    
    """
    from cobra.core import Reaction
    from copy import deepcopy
    irreversible_FVA_dict = deepcopy(irreversible_FVA_dict)
    # Popping would be easier but this creates issues with pointers in Python
    consolidated_dict = deepcopy(irreversible_FVA_dict)
    for the_reaction_id in irreversible_FVA_dict.keys():
       
        the_reaction = cobra_model.reactions.get_by_id(the_reaction_id)
        if 'reflection' in dir(the_reaction):
            if type(the_reaction.reflection) == Reaction:
                the_reflection_id = the_reaction.reflection.id
                if the_reflection_id.find("_reverse") > -1:
                    min2 = (irreversible_FVA_dict[the_reflection_id]['minimum'])
                    max2 = (irreversible_FVA_dict[the_reflection_id]['maximum']) 
                    min1 = (irreversible_FVA_dict[the_reaction_id]['minimum']) 
                    max1 = (irreversible_FVA_dict[the_reaction_id]['maximum'])
                    min2_status = (irreversible_FVA_dict[the_reflection_id]['minimum_status'])
                    max2_status = (irreversible_FVA_dict[the_reflection_id]['maximum_status']) 
                    min1_status = (irreversible_FVA_dict[the_reaction_id]['minimum_status']) 
                    max1_status = (irreversible_FVA_dict[the_reaction_id]['maximum_status'])
                    minlist = []
                    maxlist = []
                    minstatus = []
                    maxstatus = []
                    # It's possible no solutions were found when forcing the
                    # flux to be negative or positive so check to
                    # make sure numerics were found
                    if ((max1 != None) & (max1_status != 'infeasible')):
                        maxlist.append(max1)
                        maxstatus.append(max1_status)
                    if ((min1 != None) & (min1_status != 'infeasible')):
                        minlist.append(min1)
                        minstatus.append(min1_status)                        
                    if ( (min2 != None) & (min2_status != 'infeasible')):
                        maxlist.append(-1 * min2)
                        maxstatus.append(min2_status)
                    if ((max2 != None) & (max2_status != 'infeasible')):
                        minlist.append(-1 * max2)
                        minstatus.append(max2_status)
                    del consolidated_dict[the_reflection_id]
                    # also get the bounds, used to enforce consistency
                    upper_bound = max(the_reaction.upper_bound, (-1 * the_reaction.reflection.lower_bound))
                    lower_bound = min(the_reaction.lower_bound, (-1 * the_reaction.reflection.upper_bound))
                    
                    if len(minlist) > 0:
                        cur_min = min(minlist)
                        cur_index = minlist.index(min(minlist))
                        consolidated_dict[the_reaction_id]['minimum'] = cur_min
                        consolidated_dict[the_reaction_id]['minimum_status'] = minstatus[cur_index]
                    else:
                        consolidated_dict[the_reaction_id]['minimum'] = None
                        consolidated_dict[the_reaction_id]['minimum_status'] = 'infeasible'
                    if len(maxlist) > 0:
                        cur_max = max(maxlist)
                        cur_index = maxlist.index(max(maxlist))
                        consolidated_dict[the_reaction_id]['maximum'] = cur_max
                        consolidated_dict[the_reaction_id]['maximum_status'] = maxstatus[cur_index]
                    else:
                        consolidated_dict[the_reaction_id]['maximum'] = None
                        consolidated_dict[the_reaction_id]['maximum_status'] = 'infeasible'
                    # Also enforce numerical consistency
                    if enforce_consistency:
                        consolidated_dict[the_reaction_id] = check_bounds_consistency(consolidated_dict[the_reaction_id], [lower_bound, upper_bound])
    return consolidated_dict


def convert_objective_to_constraint(cobra_model, **kwargs):
    """ This function converts the objective to a turnover metabolite and
    reaction, applying the appropriate constraints.  Note if you
    wish to leave the target optimum unbounded and
    only restrict the minimum you should manipulate
    the appropriate bound.

    Arguments:
     cobra_model

    kwargs:

     objective_sense: 'maximize' or 'minimize'
      will also affect how fraction of optimum is applied

     fraction_of_optimum: fraction of the best achievable value to restrict
      the solution to.  Restricted to [0, 1].

     copy_model: whether to modify the model in place
      or return a copy.  Booelean: True/False
     
     new_reaction_name: what to call the new objective

     bound_best_optimum: a bound will always be placed at the worst
      objective value, as determined by fraction_of_optimum.
      This Boolean variable indicates whether to also add
      one to reflect the best.  True cements the best
      bound into the model and is a handy reference.
      If set to False, and the implicit bound (0 or 1000)
      returned when setting up the reaction is
      greater, then the bound will not be further constrained.

     solver parameters: solver, tolerance_feasibility,  tolerance_optimality,  tolerance_integer

    
    """
    from copy import deepcopy
    from cobra.core.Metabolite import Metabolite
    from cobra.core.Reaction import Reaction
    from numpy import array
    from math import ceil
    from math import floor

    if 'objective_sense' in kwargs:
        objective_sense = kwargs['objective_sense']
    else:
        objective_sense = 'maximize'

    if 'fraction_of_optimum' in kwargs:
        fraction_of_optimum = kwargs['fraction_of_optimum']
    else:
        fraction_of_optimum = 1.

    if 'copy_model' in kwargs:
        copy_model = kwargs['copy_model']
    else:
        copy_model = True

    if 'new_reaction_name' in kwargs:
        new_reaction_name = kwargs['new_reaction_name']
    else:
        new_reaction_name = "objective"

    if 'bound_best_optimum' in kwargs:
        bound_best_optimum = kwargs['bound_best_optimum']
    else:
        bound_best_optimum = True

    solver, tolerance_optimality, tolerance_feasibility, tolerance_barrier, tolerance_integer = get_solver_parameters(**kwargs)
    
    continue_conversion = False

    if copy_model:
        cobra_model = deepcopy(cobra_model)

    gim3e_optimize(cobra_model, objective_sense = objective_sense, solver = solver,
        tolerance_optimality = tolerance_optimality,
        tolerance_feasibility = tolerance_feasibility,
        tolerance_barrier = tolerance_barrier,
        tolerance_integer = tolerance_integer)
    wt_solution = cobra_model.solution.f

    # Just pick a new name in case we have multiple
    # Calls to this function.
    named_objective_constraint = False
    while named_objective_constraint == False:
        same_name_list = cobra_model.metabolites.query(new_reaction_name)
        if len(same_name_list) < 1:
            named_objective_constraint = True
        else:
            new_reaction_name = new_reaction_name + "_1"

    objective_metabolite = Metabolite(new_reaction_name)
    objective_reaction = Reaction(new_reaction_name)
    old_obj_list = []
    old_obj_coefficients = []
    for x in cobra_model.reactions:
        if x.objective_coefficient != 0:
            x.add_metabolites({objective_metabolite: x.objective_coefficient})
            old_obj_list.append(x)
            old_obj_coefficients.append(x.objective_coefficient)
            x.objective_coefficient = 0
            continue_conversion = True
    if continue_conversion == False:
        print("No objective detected, exiting conversion of objective to constraint...")
    else:
        objective_reaction.add_metabolites({objective_metabolite:-1})
        cobra_model.add_reactions(objective_reaction)
        if objective_sense == 'maximize':
            # Make sure the new objective isn't overly constrained
            if bound_best_optimum | (objective_reaction.upper_bound < wt_solution):
                objective_reaction.upper_bound = wt_solution            
            if wt_solution >= 0:
                objective_reaction.lower_bound = wt_solution * fraction_of_optimum
            else:
                if fraction_of_optimum > 0:
                    # See the "minimize" case for comments.
                    objective_reaction.lower_bound = wt_solution / fraction_of_optimum
                else:
                    objective_reaction.lower_bound =\
                    floor(sum(array([(old_obj_coefficients[index] * max(0, x.upper_bound))
                               for index, x in enumerate(old_obj_list) if
                               (old_obj_coefficients[index] < 0)])) +\
                    sum(array([(old_obj_coefficients[index] * min(0, x.lower_bound))
                               for index, x in enumerate(old_obj_list) if
                               (old_obj_coefficients[index] > 0)])))   
        elif objective_sense == 'minimize':
            if bound_best_optimum | (objective_reaction.lower_bound > wt_solution):
                objective_reaction.lower_bound = wt_solution            
            if wt_solution >= 0:
                # There are several options to scale the upper bound  for poor 
                # consistency in this scenario.  We could find the global worst 
                # case and scale over this range, but then the scaling becomes 
                # somewhat arbitrary.  We also could scale by 1 + (1 - 
                # fr_optimum), but then to explore ranges much worse than the 
                # optimum we need negative fr_optimum values.  A ratio is
                # therefore applied in the case fr_optimum is between 0 and 1.
                if fraction_of_optimum > 0:
                    objective_reaction.upper_bound = wt_solution / fraction_of_optimum
                # Otherwise, the biggest value is determined by the upper_bounds
                # (to avoid divide by zero errors)
                else:
                    objective_reaction.upper_bound =\
                    ceil(sum(array([(old_obj_coefficients[index] * max(0, x.upper_bound))
                               for index, x in enumerate(old_obj_list) if
                               (old_obj_coefficients[index] > 0)])) +\
                    sum(array([(old_obj_coefficients[index] * min(0, x.lower_bound))
                               for index, x in enumerate(old_obj_list) if
                               (old_obj_coefficients[index] < 0)])))    
            else:
                objective_reaction.upper_bound = wt_solution * fraction_of_optimum
        objective_reaction.objective_coefficient = 1.
        # run one optimization to store the best value
        # in cobra_model.solution
        gim3e_optimize(cobra_model, solver=solver,
            objective_sense=objective_sense,
            # Not supplying a problem and
            # have already updated the
            # model objective coefficients
            # e.g., no new_objective = {objective_reaction: 1.},
            tolerance_optimality = tolerance_optimality,
            tolerance_feasibility = tolerance_feasibility,
            tolerance_barrier = tolerance_barrier,
            tolerance_integer = tolerance_integer)
        objective_reaction.objective_coefficient = 0  

    return cobra_model


def add_sense_aware_bounds(cobra_model, **kwargs):
    """ This function adds minimal bounds to a singular objective
    to ensure future solutions will meet constraints.  Numerically,
    this can help as opposed to adding a new reaction
    to impose the bounds, which is done for
    multi-component objectives.  The model is modified in-place

	kwargs:
     objective_sense: 'maximize' or 'minimize'
      will also affect how fraction of optimum is applied

     fraction_of_optimum: fraction of the best achievable value to restrict
      the solution to

     bound_best_optimum: a bound always placed on the worst
      objective value, as determined by fraction_of_optimum.
      This Boolean variable indicates whether to also add
      one to reflect the best.

     solver parameters: solver, tolerance_feasibility,  tolerance_optimality,  tolerance_integer
	 
	returns:
	 cobra_model if effective, None if not
	 

    
    """
    import types
    from copy import deepcopy
    from cobra.core.Metabolite import Metabolite
    from cobra.core.Reaction import Reaction
    from numpy import array
    from math import ceil
    from math import floor

    if 'objective_sense' in kwargs:
        objective_sense = kwargs['objective_sense']
    else:
        objective_sense = 'maximize'

    if 'fraction_of_optimum' in kwargs:
        fraction_of_optimum = kwargs['fraction_of_optimum']
    else:
        fraction_of_optimum = 1.

    if 'bound_best_optimum' in kwargs:
        bound_best_optimum = kwargs['bound_best_optimum']
    else:
        bound_best_optimum = True

    solver, tolerance_optimality, tolerance_feasibility, tolerance_barrier, tolerance_integer = get_solver_parameters(**kwargs)

    continue_conversion = True

    nonzero_objective_list = [x.id for x in cobra_model.reactions if abs(x.objective_coefficient) > 0]

    if len(nonzero_objective_list) != 1:
        print("Only call add_sense_aware_bounds for objectives with exactly one component, exiting ...")
    else:
        # Have noted overconstraining solution may present
        # challenges in finding the optimum.  This is not
        # mathematically needed but is numerically helpful
        gim3e_optimize(cobra_model, objective_sense = objective_sense, solver = solver,
            tolerance_optimality = tolerance_optimality,
            tolerance_feasibility = tolerance_feasibility,
            tolerance_barrier = tolerance_barrier,
            tolerance_integer = tolerance_integer)
        
        wt_solution = cobra_model.solution.f
        objective_reaction = cobra_model.reactions.get_by_id(nonzero_objective_list[0])
        # Usually the coefficient will be 1.
        if type(wt_solution) != types.NoneType:
			best_obj_reaction_value = wt_solution / objective_reaction.objective_coefficient
        else:
			continue_conversion = False
			print("Unable to optimize input model while adding sense aware bounds.")
			return None
	
        if (objective_sense == 'maximize') & continue_conversion:
            # Since we aren't making a new objective function,
            # we respect existing limits on the objective
            # and don't include a check on
            # objective_reaction.upper_bound < best_obj_reaction_value
            if bound_best_optimum:
                objective_reaction.upper_bound = best_obj_reaction_value            
            if best_obj_reaction_value >= 0:
                if fraction_of_optimum > 0:
                    proposed_lower_bound = best_obj_reaction_value * fraction_of_optimum
                    if objective_reaction.lower_bound < proposed_lower_bound:
                        objective_reaction.lower_bound = proposed_lower_bound
            else:
                if fraction_of_optimum > 0:
                    proposed_lower_bound = best_obj_reaction_value / fraction_of_optimum
                    if objective_reaction.lower_bound < proposed_lower_bound:
                        # See the "minimize" case for comments.
                        objective_reaction.lower_bound = proposed_lower_bound
                # else: if we are going to maximize & best_value is negative,
                # if the fr of optimum is 0 we don't need to add a bound.

        elif (objective_sense == 'minimize') & (continue_conversion):
            if bound_best_optimum:
            # Since we aren't making a new objective function,
            # we respect existing limits on the objective            
				objective_reaction.lower_bound = best_obj_reaction_value            
            if best_obj_reaction_value >= 0:
                # There are several options to scale the upper bound  for poor 
                # consistency in this scenario.  We could find the global worst 
                # case and scale over this range, but then the scaling becomes 
                # somewhat arbitrary.  We also could scale by 1 + (1 - 
                # fr_optimum), but then to explore ranges much worse than the 
                # optimum we need negative fr_optimum values.  A ratio is
                # therefore applied in the case fr_optimum is between 0 and 1.
                if fraction_of_optimum > 0:
                    proposed_upper_bound = best_obj_reaction_value / fraction_of_optimum
                    if objective_reaction.upper_bound > proposed_upper_bound:
                        objective_reaction.upper_bound = proposed_upper_bound
                
                # else: if fraction_of_optimum == 0 & we want to minimize,
                # then we don't
                # care how big the objective gets.
                # Then we don't adjust upper_bound
            else:
                if fraction_of_optimum > 0:
                    proposed_upper_bound = best_obj_reaction_value * fraction_of_optimum
                    if objective_reaction.upper_bound > proposed_upper_bound:
                        objective_reaction.upper_bound = proposed_upper_bound
        # run one optimization to store the best value
        # in cobra_model.solution
        gim3e_optimize(cobra_model, solver=solver,
            objective_sense=objective_sense,
            # Don't need to call since
            # objective has not changed
            # new_objective=objective_reaction,
            tolerance_optimality=tolerance_optimality,
            tolerance_feasibility=tolerance_feasibility,
            tolerance_barrier=tolerance_barrier,
            tolerance_integer = tolerance_integer)
        # Don't reset the objective coefficient to 0      

    return cobra_model


def convert_FVA_to_constraint(cobra_model, FVA_dict, **kwargs):
    """ Take ranges from FVA and converts them into model constraints.
    Potentially useful for sampling.  Best to first call 'turn_off_unused_reactions'

    Also note that this script has not been  tested
    with models not converted to irreversible, eg with negative
    flux ranges.

    Arguments:
     cobra_model
     FVA_dict: results from irreversible_flux_variability_analysis
     copy_model: if False, the model will be modified in place.

    Keyword arguments:
     enforce_original_bounds: if true, we will not set bounds outside the
      constraints in the cobra model object
     epsilon: a small number, some reactions may require very low flux
      for growth.  Should be based on solver tolerance in the FVA
      step.  1x epsilon will be used to help set the boundaries.

    Returns:
     cobra_model: model with updated bounds


    """
    from math import floor
    from math import ceil
    from copy import deepcopy
    from cobra.core.Reaction import Reaction

    if 'copy_model' in kwargs:
        copy_model = kwargs['copy_model']
    else:
        copy_model = True

    if 'enforce_original_bounds' in kwargs:
        enforce_original_bounds = kwargs['enforce_original_bounds']
    else:
        enforce_original_bounds = True

    if 'epsilon' in kwargs:
        epsilon = kwargs['epsilon']
    else:
        epsilon = 1E-7

    if copy_model:
        cobra_model = deepcopy(cobra_model)

    the_model_reaction_ids = [x.id for x in cobra_model.reactions]
    for reaction_id in FVA_dict.keys():
        if reaction_id in the_model_reaction_ids:
            the_reaction = cobra_model.reactions.get_by_id(reaction_id)
            the_max = FVA_dict[reaction_id]['maximum']
            the_min = FVA_dict[reaction_id]['minimum']
            the_max_ref = None
            the_min_ref = None
            has_ref = False

            if 'reflection' in dir(the_reaction):
                if type(the_reaction.reflection) == Reaction:
                    reflection_id = the_reaction.reflection.id
                    the_max_ref = FVA_dict[reflection_id]['maximum']
                    the_min_ref = FVA_dict[reflection_id]['minimum']
                    has_ref = True

            # We only try to constrain if at least one
            # of the forward or reverse reactions is not "none",
            # Otherwise this reaction was not treated well by
            # the solver, and we leave it alone.
            if ((the_max != None) | (the_min != None) | (the_max_ref != None) | (the_min_ref != None)):
                if (the_max != None):
                    the_max = the_max + epsilon
                    the_max = ceil(the_max / epsilon) * epsilon
                    if enforce_original_bounds == True:
                        if ((the_max < the_reaction.upper_bound) & (the_max > the_reaction.lower_bound)):
                            the_reaction.upper_bound = the_max
                    else:
                        the_reaction.upper_bound = the_max
                elif (the_min == None) & (the_max_ref != None) & (the_min_ref != None):
                    # None shows up in FVA results for irreversible
                    # models when flux is allowed only in the opposite
                    # direction or a small flux below the
                    # solver tolerance is required.  Since
                    # we know at least one of the
                    # solutions worked by the preceeding or statement,
                    # this must be a case where flux in the opposite
                    # direction only is allowed and it
                    # is safe to set the bounds equal to zero.
                    the_reaction.upper_bound = 0
    
                the_min = FVA_dict[reaction_id]['minimum']
                if (the_min != None):
                    the_min = the_min - epsilon
                    the_min = floor(the_min / epsilon) * epsilon
                    if enforce_original_bounds == True:
                        if ((the_min < the_reaction.upper_bound) & (the_min > the_reaction.lower_bound)):
                            the_reaction.lower_bound = the_min
                    else:
                        the_reaction.lower_bound = the_min                        
                elif (the_max == None) & (the_max_ref != None) & (the_min_ref != None):
                    the_reaction.lower_bound = 0                        
        else:
            print('Did not detect id '+ reaction_id +' in model.  Skipping this reaction.')
    return cobra_model


def turn_off_unused_reactions(cobra_model, **kwargs):

    """ Take the results from FVA to identify candidates where bounds on unused
    reactions in the optimal network states are set to zero.  As a check, the
    model is tested to verify the objective score is not adversely affected
    (fraction of optimum).

    Arguments:
     cobra_model: an objective (typically growth) should be specified in the model.
      The model should be irreversible.

     kwargs:
      FVA_dict: results of flux variability analysis (DO NOT CONVERT TO REVERSIBLE).
       Used to identify reactions
       with low flux for removal.  If not provided and no specific candidates are
       provided, FVA will be performed for all reactions to assist candidate
       selection.
     
      epsilon: a threshold, reactions with UB less than this value are considered
       for turning off.

      solver & tolerance parameters: solver, tolerance_optimality, tolerance_feasibility, tolerance_barrier, tolerance_integer

      candidate_reaction_id_list: optional, a list of candidate reactions to
       consider for inactivation.  If empty all low-flux reactions will be
       considered.
     
      fraction_of_optimum: a threshold for the objective reaction for the
       returned model relative to the input.
       Note only maximization problems are currently supported 

      constrained_reaction_list: optional, a list of reactions that should
       be constrained not to degrade at the optimal objective.  e.g. we can limit
       to a fraction of the best growth function value, then further constrain
       so that the removal of additional reactions does not change the range in
       penalty.

      verbose: if true, report the reactions as they are tested to track progress.

     Returns:
      cobra_model: a model with all of the ~zero flux reactions
       whose removal does not degrade the specified objective deactivated
      turned_off_reaction_list: list of reaction ids that were deactivated

      Note: currently only set up for objectives that are maximization problems.

     
    """
    from copy import deepcopy
    from cobra.core.Reaction import Reaction
    from numpy import nan

    if 'FVA_dict' in kwargs:
        FVA_dict = kwargs['FVA_dict']
    else:
        FVA_dict = {}

    if 'epsilon' in kwargs:
        epsilon = kwargs['epsilon']
    else:
        epsilon = 1E-7

    solver, tolerance_optimality, tolerance_feasibility, tolerance_barrier, tolerance_integer = get_solver_parameters(**kwargs)

    if 'candidate_reaction_id_list' in kwargs:
        candidate_reaction_id_list = kwargs['candidate_reaction_id_list']
    else:
        candidate_reaction_id_list = []

    if 'fraction_of_optimum' in kwargs:
        fraction_of_optimum = kwargs['fraction_of_optimum']
    else:
        fraction_of_optimum = 1.

    if 'constrained_reaction_list' in kwargs:
       constrained_reaction_list = kwargs['constrained_reaction_list']
    else:
        constrained_reaction_list = []

    if 'verbose' in kwargs:
        verbose = kwargs['verbose']
    else:
        verbose = False        
    
    result_dict = {}
    # Copy the model since we may modify it.
    cobra_model = deepcopy(cobra_model)
    turned_off_reaction_list = []    
    turned_off_objective_values = []
    # We need to use a tolerance <= epsilon
    # in order to properly evaluate whether
    # the reaction can be inactivated.
    if len(candidate_reaction_id_list) > 0:
        # We enforce the reversibility here, can't propose a reaction
        # for inactivation without its partner if it is reversible.
        new_FVA_dict = {}
        model_reaction_list = [x.id for x in cobra_model.reactions]
        for reaction_id in candidate_reaction_id_list:
            if (((reaction_id + "_reverse") in model_reaction_list) &
                (not(reaction_id+"_reverse" in candidate_reaction_id_list))):
                candidate_reaction_id_list.append(reaction_id+"_reverse")
            elif (((reaction_id.rstrip("_reverse")) in model_reaction_list) &
                  (not(reaction_id.rstrip("_reverse") in candidate_reaction_id_list))):
                candidate_reaction_id_list.append(reaction_id.rstrip("_reverse"))                 

    # Run FVA to ID inactivation candidates if an FVA dict was not supplied
    # and the candidates were not specified
    if (len(FVA_dict.keys()) < 1):
        if (len(candidate_reaction_id_list) > 0):
            the_FVA_reaction_ids = candidate_reaction_id_list
        else:
            the_FVA_reaction_ids = None
        new_FVA_dict = irreversible_flux_variability_analysis(cobra_model,
            fraction_of_optimum = fraction_of_optimum,
            objective_sense='maximize',
            the_reactions= the_FVA_reaction_ids,
            solver = solver,
            tolerance_optimality = tolerance_optimality,
            tolerance_feasibility = tolerance_feasibility,
            tolerance_barrier = tolerance_barrier,
            tolerance_integer = tolerance_integer)
    elif (len(candidate_reaction_id_list) > 0):
        # We ignore candidates without FVA results
        new_FVA_dict = {}
        # Potential artifact if FVA dict has already been converted to reversible.
        for candidate_reaction in candidate_reaction_id_list:
            if FVA_dict.has_key(candidate_reaction):
                new_FVA_dict[candidate_reaction] = deepcopy(FVA_dict[candidate_reaction])
    else:
        # If no candidate list is given we just use the whole set of results
        new_FVA_dict = deepcopy(FVA_dict)
    # Pick out candidate reactions for elimination based on the FVA results
    # Need to convert to reversible, also enforce numerical consistency
    if verbose:
        print("Constructed an irreversible FVA dict with " + str(len(new_FVA_dict.keys())) + " irreversible reaction keys to query for inactivation.")    

    FVA_reversible_dict = convert_FVA_dict_to_reversible(new_FVA_dict, cobra_model)
    if verbose:
        print("Constructed a reversible FVA dict with " + str(len(FVA_reversible_dict.keys())) + " reaction keys to query for inactivation.")

    elimination_candidates = []
    for reaction_id in FVA_reversible_dict.keys():
        # we consider a reaction for inactivation if none of its
        # fluxes is larger than epsilon.  Leave infeasible
        # alone.
        magnitudes = [0]
        if ((FVA_reversible_dict[reaction_id]['maximum_status'] != 'infeasible')
            & (FVA_reversible_dict[reaction_id]['maximum'] != None)):
            magnitudes.append(abs(FVA_reversible_dict[reaction_id]['maximum']))
        if ((FVA_reversible_dict[reaction_id]['minimum_status'] != 'infeasible')
            & (FVA_reversible_dict[reaction_id]['minimum'] != None)):
            magnitudes.append(abs(FVA_reversible_dict[reaction_id]['minimum']))
        if max(magnitudes) <= epsilon:
            set_bounds_magnitudes = []
            the_reaction = cobra_model.reactions.get_by_id(reaction_id)
            set_bounds_magnitudes.extend([abs(the_reaction.upper_bound),
                                          abs(the_reaction.lower_bound)])
            if 'reflection' in dir(the_reaction):
                if type(the_reaction.reflection) == Reaction:
                    set_bounds_magnitudes.extend([abs(the_reaction.reflection.upper_bound),
                                                  abs(the_reaction.reflection.lower_bound)])
            # Ignore if the bounds are already zero
            if max(set_bounds_magnitudes) > 0:
                elimination_candidates.append(reaction_id)

    if verbose:
        print("Identified " + str(len(elimination_candidates)) + " reactions as candidates for deactivation based on FVA.")

    # If specified, we need to establish that additional reaction ranges do not
    # leave the original bounds at the optimum following the elimination of a
    # reaction, so we incorporate these as additional constraints.
    # This ensures the range of new values after removing the reaction fall
    # within the old allowed range.  This leaves the possibility that there
    # may be some degradation of, for example, the best achievable total penalty.
    # To be fully conservative - e.g. if there are competing
    # biological interests, we can re-run this script multiple times, switching
    # the constraints and the objective and use the intersection of inactivation
    # candidates for inactivation
    constraint_FVA_dict = {}
    original_boundary_dict = {}                    
    if (len(constrained_reaction_list) > 0):
        # First check to make sure forward and reverse reaction pairs are included
        model_reaction_list = [x.id for x in cobra_model.reactions]
        for reaction_id in constrained_reaction_list:
            if (((reaction_id + "_reverse") in model_reaction_list) &
                (not(reaction_id+"_reverse" in constrained_reaction_list))):
                constrained_reaction_list.append(reaction_id+"_reverse")
            elif (((reaction_id.rstrip("_reverse")) in model_reaction_list) &
                  (not(reaction_id.rstrip("_reverse") in constrained_reaction_list))):
                constrained_reaction_list.append(reaction_id.rstrip("_reverse"))                 

        constraint_FVA_dict = irreversible_flux_variability_analysis(cobra_model,
                              fraction_of_optimum = fraction_of_optimum,
                              objective_sense='maximize',
                              the_reactions = constrained_reaction_list,
                              solver=solver,
                              tolerance_optimality = solver_tolerance,
                              tolerance_feasibility = solver_tolerance,
                              tolerance_barrier = 0.0001*solver_tolerance,
                              tolerance_integer = 1E-9,
                              error_reporting = None,
                              number_of_processes = 1)
        for constrained_reaction_id in constrained_reaction_list:
            # if the constrained reactions truly are very constrained then
            # we need to relax the result somewhat so the solver will
            # function
            the_reaction = cobra_model.reactions.get_by_id(constrained_reaction_id)
            maximum_status = constraint_FVA_dict[constrained_reaction_id]['maximum_status']
            minimum_status = constraint_FVA_dict[constrained_reaction_id]['minimum_status']
            if ((maximum_status in acceptable_solution_strings) & (minimum_status in acceptable_solution_strings)):
                current_maximum = constraint_FVA_dict[constrained_reaction_id]['maximum']
                current_minimum = constraint_FVA_dict[constrained_reaction_id]['minimum']
                if (abs(current_maximum - current_minimum) < solver_tolerance):
                    constraint_FVA_dict[constrained_reaction_id]['maximum'] = max(
                        min((current_maximum + solver_tolerance),
                            the_reaction.upper_bound),the_reaction.lower_bound)
                    constraint_FVA_dict[constrained_reaction_id]['minimum'] = min(
                        max((current_minimum - solver_tolerance),
                            the_reaction.lower_bound),the_reaction.upper_bound)
            original_boundary_dict[constrained_reaction_id] = {}
            original_boundary_dict[constrained_reaction_id]['upper_bound'] = deepcopy(
                the_reaction.upper_bound)
            original_boundary_dict[constrained_reaction_id]['lower_bound'] = deepcopy(
                the_reaction.lower_bound)            
            if maximum_status in acceptable_solution_strings:
                if constraint_FVA_dict[constrained_reaction_id]['maximum'] < the_reaction.upper_bound:
                    the_reaction.upper_bound = deepcopy(constraint_FVA_dict[constrained_reaction_id]['maximum'])
            if minimum_status in acceptable_solution_strings:
                if constraint_FVA_dict[constrained_reaction_id]['minimum'] > the_reaction.lower_bound:
                    the_reaction.lower_bound = deepcopy(constraint_FVA_dict[constrained_reaction_id]['minimum'])

    # Solve one more time with the new constraint imposed
    gim3e_optimize(cobra_model, objective_sense = 'maximize', 
                   the_problem = None, solver = solver,  
                   error_reporting = None,
                   tolerance_optimality = tolerance_optimality, 
                   tolerance_feasibility = tolerance_feasibility,
                   tolerance_barrier = tolerance_barrier,
                   tolerance_integer = tolerance_integer)
    
    objective_components = [x.id for x in cobra_model.reactions if (x.objective_coefficient != 0)]
    if len(objective_components) > 1:
        print("Warning, multi-component objectives not yet supported in 'turn_off_unused_reactions()'.")
    else:
        objective_reaction_id = objective_components[0]

    # Now we "bake in" the optimization constraint, this is
    # potentially more rigorous than incorporating a >= type of comparison
    # afterwards and we can drop the allowance for a comparison tolerance
    # since the solver handles it implicitly.
    the_objective_reaction = cobra_model.reactions.get_by_id(objective_reaction_id)
    initial_objective_reaction_lower_bound = the_objective_reaction.lower_bound
    the_objective_reaction.lower_bound = fraction_of_optimum * cobra_model.solution.f
    initial_best_solution = gim3e_optimize(cobra_model, objective_sense = 'maximize', 
                   the_problem = None, solver = solver,  
                   error_reporting = None,
                   tolerance_optimality = tolerance_optimality, 
                   tolerance_feasibility = tolerance_feasibility,
                   tolerance_barrier = tolerance_barrier,
                   tolerance_integer = tolerance_integer)
    best_objective = cobra_model.solution.f

    # Back up the cobra_model so we have an original to reference
    # reaction bounds, etc...
    cobra_model_initial_opt = deepcopy(cobra_model)    

    # Verify the reactions can be turned off, one at a time
    reaction_counter = 1
    for reaction_id in elimination_candidates:
        if verbose:
            print("Testing reaction "+ str(reaction_counter) + " of " +
              str(len(elimination_candidates)) + " for inactivation.")
        reaction_counter = reaction_counter + 1
        the_solution = initial_best_solution
        the_reaction = cobra_model.reactions.get_by_id(reaction_id)
        original_upper_bound = the_reaction.upper_bound
        original_lower_bound = the_reaction.lower_bound
        the_reaction.upper_bound = 0
        the_reaction.lower_bound = 0        
        the_index = cobra_model_initial_opt.reactions.index(
            cobra_model_initial_opt.reactions.get_by_id(the_reaction.id))
        if (cobra_model_initial_opt.solution.x[the_index] != 0.):
            the_solution = None          
        has_ref = False
        if 'reflection' in dir(the_reaction):
            if type(the_reaction.reflection) == Reaction:
                has_ref = True
                the_reflection_id = the_reaction.reflection.id
                original_ref_upper_bound = the_reaction.reflection.upper_bound
                original_ref_lower_bound = the_reaction.reflection.lower_bound
                the_reaction.reflection.upper_bound = 0
                the_reaction.reflection.lower_bound = 0
                the_index = cobra_model_initial_opt.reactions.index(
                    cobra_model_initial_opt.reactions.get_by_id(the_reflection_id))
                if (cobra_model_initial_opt.solution.x[the_index] != 0.):
                    the_solution = None
        gim3e_optimize(cobra_model, objective_sense = 'maximize', 
                   the_problem = None, solver = solver,  
                   error_reporting = None,
                   tolerance_optimality = tolerance_optimality, 
                   tolerance_feasibility = tolerance_feasibility,
                   tolerance_barrier = tolerance_barrier,
                   tolerance_integer = tolerance_integer)
        # Revert the bounds
        the_reaction.upper_bound = original_upper_bound
        the_reaction.lower_bound = original_lower_bound
        if has_ref:
            the_reaction.reflection.upper_bound = original_ref_upper_bound
            the_reaction.reflection.lower_bound = original_ref_lower_bound  
        if cobra_model.solution.status in acceptable_solution_strings:
            test_obj = cobra_model.solution.f
            # We incorporated the limits as a constraint so if it solved,
            # it is OK to test this reaction for being turned off.
            turned_off_objective_values.append(test_obj)
            turned_off_reaction_list.append(reaction_id)            

    # Now we inactivate all candidates simultaneously
    # Note we have these backed up in cobra_model_initial_opt
    for the_reaction_id in turned_off_reaction_list:
        the_reaction = cobra_model.reactions.get_by_id(the_reaction_id)
        the_reaction.upper_bound = 0
        the_reaction.lower_bound = 0        
        if 'reflection' in dir(the_reaction):
            if type(the_reaction.reflection) == Reaction:
                has_ref = True
                the_reaction.reflection.upper_bound = 0
                the_reaction.reflection.lower_bound = 0

    # Now that the candidates have been identified we loop through again until
    # we can find identify a larger set of eliminated reactions that 
    # best maintain growth (and constrained reactions together).  The search heuristic is
    # simply to reactivate reactions that have the largest impact on the full model
    # individually.
    
    reactivated_counter = 0
    the_solution_status = 'infeasible'
    while the_solution_status not in acceptable_solution_strings:
        gim3e_optimize(cobra_model, objective_sense = 'maximize', 
                   the_problem = None, solver = solver,  
                   error_reporting = None,
                   tolerance_optimality = tolerance_optimality, 
                   tolerance_feasibility = tolerance_feasibility,
                   tolerance_barrier = tolerance_barrier,
                   tolerance_integer = tolerance_integer)      
        if cobra_model.solution.status not in acceptable_solution_strings:
            if len(turned_off_objective_values) > 0:
                reactivated_counter = reactivated_counter + 1
                next_to_go_index = turned_off_objective_values.index(min(turned_off_objective_values))
                turned_off_objective_values.pop(next_to_go_index)
                the_reaction_id = turned_off_reaction_list.pop(next_to_go_index)
                the_reaction = cobra_model.reactions.get_by_id(the_reaction_id)
                the_reaction.upper_bound = cobra_model_initial_opt.reactions.get_by_id(the_reaction.id).upper_bound
                the_reaction.lower_bound = cobra_model_initial_opt.reactions.get_by_id(the_reaction.id).lower_bound
                if 'reflection' in dir(the_reaction):
                    if type(the_reaction.reflection) == Reaction:
                        the_reaction.reflection.upper_bound = cobra_model_initial_opt.reactions.get_by_id(the_reaction.id).reflection.upper_bound
                        the_reaction.reflection.lower_bound = cobra_model_initial_opt.reactions.get_by_id(the_reaction.id).reflection.lower_bound
                if verbose:
                    print("Could not find acceptable solution, reactivating reaction " + str(reactivated_counter) + ".")
            else:
                if verbose:
                    print("Could not find an acceptable solution when removing reactions, all reactions are reactivated.")
                the_solution_status = 'optimal'
        else:
            the_solution_status = 'optimal'
        

    # We may need to revert reaction bounds
    if (len(constrained_reaction_list) > 0):
        for constrained_reaction_id in constrained_reaction_list:
            the_reaction = cobra_model.reactions.get_by_id(constrained_reaction_id)            
            the_reaction.upper_bound = deepcopy(original_boundary_dict[constrained_reaction_id]['upper_bound'])
            the_reaction.lower_bound = deepcopy(original_boundary_dict[constrained_reaction_id]['lower_bound'])

    cobra_model.reactions.get_by_id(objective_reaction_id).lower_bound = initial_objective_reaction_lower_bound
    
    return (cobra_model, turned_off_reaction_list)


def check_bounds_consistency(tmp_dict, the_reaction):
    """ A simple function to ensure the FVA bounds do not conflict.  Widens the
    range a bit to accomodate the solver tolerance if there is a conflict.
    
    Arguments:
     tmp_dict: dictionary with the minimum, maximum, solver status associated
      with the minimum, and solver status for the maximum.
     the_reaction: either a reaction object or list of bounds
     
    Output:
     tmp_dict: input modified in place and also returned


    """
    
    if type(the_reaction) == list:
        lower_bound = min(the_reaction)
        upper_bound = max(the_reaction)
    else:
        lower_bound = the_reaction.lower_bound
        upper_bound = the_reaction.upper_bound

    the_min = tmp_dict['minimum']
    the_max = tmp_dict['maximum']
    the_min_status = tmp_dict['minimum_status']
    the_max_status = tmp_dict['maximum_status']

    # First, the min can't exceed the max, this can happen
    # with FVA results near zero
    # Rather than loosening the range we simply recorrect
    # the ordering.  Tolerances are added in e.g.
    # when resetting model bounds.  It is better
    # to maintain the actual solver values here.
    if ((the_min != None) & (the_max != None)):
        if the_min > the_max:
            tmp_dict['minimum'] = the_max
            tmp_dict['maximum'] = the_min
            tmp_dict['minimum_status'] = the_max_status
            tmp_dict['maximum_status'] = the_min_status
    # Preceeding might fail and might need to revert to adding in a tolerance

    # Second, the solver may report '-0'
    # Numerically valid according to Python, but this is not desirable
    for the_key in ['minimum', 'maximum']:
        if tmp_dict[the_key] != None:
            if tmp_dict[the_key] == 0:
                tmp_dict[the_key] = 0.0
        

    # Third neither the max or min can be outside the bounds if they are reported as.
    # optimal solutions.  This is just "trimming" the solutions to be compatible
    # with downstream processing
    if (tmp_dict['maximum_status'] in acceptable_solution_strings):
        if (tmp_dict['maximum'] > upper_bound):
            # Took out earlier check to verify id this is above but not within
            # tolerance... not a valid scenario.
            tmp_dict['maximum'] = upper_bound
        elif (tmp_dict['maximum'] < lower_bound):
            tmp_dict['maximum'] = lower_bound
    if (tmp_dict['minimum_status'] in acceptable_solution_strings):
        if (tmp_dict['minimum'] < lower_bound):
            tmp_dict['minimum'] = lower_bound
        elif (tmp_dict['minimum'] > upper_bound):
            tmp_dict['minimum'] = upper_bound
            
    return tmp_dict


def gim3e_optimize(cobra_model, solver='cplex', error_reporting=True, **kwargs):
    """ A variation on optimize, constructed for conservative numerical
    calculations.

    Arguments:
     cobra_model: the COBRA model object to import
     solver: preferred solver
     error_reporting: whether to ouput optimization
      errors (e.g. due to infeasibilities, etc...)
     kwargs: optional keyword arguments.  Note that
      because glpk requires its tolerances
      whenever it is run, updates are applied to all
      'the_solution'/lp objects for consistency
      for all solvers is performed

    """
    from cobra.flux_analysis.objective import update_objective
    from cobra import solvers
    from copy import deepcopy
    import types

    ### Solver specific parameters
    from cobra.solvers.parameters import status_dict, \
        parameter_mappings, parameter_defaults

    solver = check_solver(solver)

    # For transparency, declare parameter defaults locally
    # These are overwitten by any kwargs as configuration_parameters
    parameter_defaults = deepcopy(parameter_defaults[solver])
    parameter_defaults = {'new_objective': None,
                          'objective_sense': 'maximize',
                          'min_norm': 0,
                          'the_problem': None, 
                          'tolerance_optimality': 1e-7,
                          'tolerance_feasibility': 1e-7,
                          'tolerance_integer': 1e-9, 
                          'error_reporting': None,
                          'print_solver_time': False,
                          'quadratic_component': None}
    
    # Also, pop integer gap defaults if any
    if 'MIP_gap_abs' in parameter_defaults.keys():
        parameter_defaults.pop('MIP_gap_abs')
    elif 'MIP_gap' in parameter_defaults.keys():
        parameter_defaults.pop('MIP_gap')

    # Separate out parameters that should only be set for mips,
    # that give cplex problems if set for lps
    mip_only_parameter_defaults = {'tolerance_barrier': 1E-11}
    mip_only_parameter_dict = {}
    for the_parameter in mip_only_parameter_defaults.keys():
        if the_parameter in kwargs.keys():
            mip_only_parameter_dict[the_parameter] = kwargs[the_parameter]
            kwargs.pop(the_parameter)
        else:
            mip_only_parameter_dict[the_parameter] = mip_only_parameter_defaults[the_parameter]

    # Assume we have addressed potential solver issues in the calling function
    the_solver = solvers.solver_dict[solver]
    if 'the_problem' not in kwargs:
        kwargs['the_problem'] = None

    # Update objectives if they are new.
    # update objective will be looking reaction objects
    # or reaction integer indices as reaction objects or integers
    # Convert id's
    reaction_id_list = []
    if 'new_objective' in kwargs and \
           kwargs['new_objective'] not in ['update problem', None]:
        new_objectives = {}
        if isinstance(kwargs['new_objective'], str):
            reaction_id_list = [x.id for x in cobra_model.reactions]
            if kwargs['new_objective'] in reaction_id_list:
                the_key = cobra_model.reactions.get_by_id(kwargs['new_objective'])
                new_objectives[the_key] = 1
                kwargs['new_objective'] = new_objectives
        elif isinstance(kwargs['new_objective'], dict):
            for the_objective in kwargs['new_objective'].keys():
                if isinstance(the_objective, str):
                    if len(reaction_id_list) == 0:
                        reaction_id_list = [x.id for x in cobra_model.reactions]
                    if the_objective in reaction_id_list:
                        the_key = cobra_model.reactions.get_by_id(the_objective)
                        new_objectives[the_key] = kwargs['new_objective'][the_objective]
                else:
                    new_objectives[the_objective] = kwargs['new_objective'][the_objective]
            kwargs['new_objective'] = new_objectives
        elif isinstance(kwargs['new_objective'], list):
            for the_objective in kwargs['new_objective']:
                if isinstance(the_objective, str):
                    if len(reaction_id_list) == 0:
                        reaction_id_list = [x.id for x in cobra_model.reactions]
                    if the_objective in reaction_id_list:
                        the_key = cobra_model.reactions.get_by_id(the_objective)
                        new_objectives[the_key] = 1.
                else:
                    new_objectives[the_objective] = 1.
            kwargs['new_objective'] = new_objectives     
        update_objective(cobra_model, kwargs['new_objective'])    

    alt_cplex_flag = False
    alt_gurobi_flag = False    

    # If not given, construct the lp items manually so we can
    # supply custom settings.
    # Also note that if supplied, the_problem gets updated here
    # Note the call to update the problem just adjusts
    # the objective and bounds; any additional parameter
    # alterations are made on the call to solve_problem
    if solver == 'cplex':
        the_methods = [1, 2, 3, 4, 5, 6]
        parameter_defaults.update({'lp_method': 1,
                    'tolerance_markowitz': 0.9})
                    # It can be detrimental to set the gap too small
                    # 'MIP_gap_abs': 0,
                    # 'MIP_gap': 0})
        configuration_parameters = deepcopy(parameter_defaults)
        configuration_parameters.update(kwargs)
        # Allow a little bit of slack in the integer solution
        # based on the selected tolerance
        # As a rule of thumb, sall absolute gaps are pretty hard for the solver
        # But seem to be able to find good solutions with 1E-6.
        #configuration_parameters['MIP_gap_abs'] = max(1E-6, configuration_parameters['tolerance_feasibility'])
        #configuration_parameters['MIP_gap'] = 0        
        if configuration_parameters['the_problem'] == None:
            the_problem = the_solver.create_problem(cobra_model, **configuration_parameters)
            the_problem.parameters.read.scale.set(-1)
            the_problem.parameters.emphasis.numerical.set(1)
            the_problem.parameters.preprocessing.presolve.set(0)
            alt_cplex_flag = True
        else:
            alt_cplex_flag = True            
            the_problem = configuration_parameters['the_problem']
            the_solver.update_problem(the_problem, cobra_model, **configuration_parameters)

    elif solver == 'glpk':
        the_methods = [1, 2, 3]
        configuration_parameters = deepcopy(parameter_defaults)
        configuration_parameters.update(kwargs)        
        # Might want to adjust tm_lim        
        if configuration_parameters['the_problem'] == None:
            # Note most GLPK tolerance parameters are set at execution
            # and not lp construction
            the_problem = the_solver.create_problem(cobra_model, **configuration_parameters)
            # Additional arguments here if needed
        else:
            the_problem = configuration_parameters['the_problem']
            the_solver.update_problem(the_problem, cobra_model, **configuration_parameters)

    # Gurobi and java glpk were not tested extensively.
    # Add your own tweaks here.
    elif solver == 'gurobi':      
        the_methods = [0, 2, 1]
        parameter_defaults.update({'lp_method': 1,
                    'tolerance_markowitz': 1E-4})
                    # 'MIP_gap_abs': 0,
                    # 'MIP_gap': 0})
        configuration_parameters = deepcopy(parameter_defaults)
        configuration_parameters.update(kwargs)
        # Allow a little bit of slack in the integer solution
        # based on the selected tolerance        
        configuration_parameters['MIP_gap_abs'] = configuration_parameters['tolerance_feasibility']
        configuration_parameters['MIP_gap'] = 0 
        if configuration_parameters['the_problem'] == None:
            the_problem = the_solver.create_problem(cobra_model, **configuration_parameters)
            the_problem.setParam('presolve', 0)
            the_problem.setParam('scaleflag', 0)
            # the_problem.setParam('mipgapabs', 0)
            alt_gurobi_flag = True            
            # Additional arguments here if needed
        else:
            the_problem = configuration_parameters['the_problem']
            alt_gurobi_flag = True            
            the_solver.update_problem(the_problem, cobra_model, **configuration_parameters)

    # didn't see this parameter in the cplex or gurobi dict
    if solver == 'cplex':
        the_problem.parameters.mip.tolerances.integrality.set(configuration_parameters['tolerance_integer'])
    if solver == 'gurobi':
        the_problem.setParam('intfeastol', configuration_parameters['tolerance_integer'])

    # Only try to set troublesome parameters for a MILP
    # in cplex or it will return an error message
    if solver == 'cplex':
        problem_type = the_problem.problem_type[the_problem.get_problem_type()]
        if problem_type != 'LP':
            the_solver.update_problem(the_problem, cobra_model, **mip_only_parameter_dict)
    else:
        the_solver.update_problem(the_problem, cobra_model, **mip_only_parameter_dict)

    ###Try to solve the problem using other methods if the first method doesn't work
    lp = the_problem
    try:
        lp_method = kwargs['lp_method']
    except:
        lp_method = 1
    if lp_method in the_methods:
        the_methods.remove(lp_method)
    #Start with the user specified method
    the_methods.insert(0, lp_method)
    # lp.set_log_stream("cplex_log.txt", lambda a:  a + " LOG " )
    # lp.set_warning_stream("cplex_warnings.txt", lambda a:  " WARNING " + a)
    # lp.set_error_stream("cplex_errors.txt", lambda a:  " ERROR " + a)
    # lp.set_results_stream("cplex_results.txt", lambda a:  " RESULT " + a)    
    acceptable_solution = None
    acceptable_solution_status = None
    for the_method in the_methods:
        configuration_parameters['lp_method'] = the_method
        try:
            if not (alt_cplex_flag | alt_gurobi_flag):
                status = the_solver.solve_problem(lp, **configuration_parameters)
            elif alt_cplex_flag:
                status = local_solve_cplex_problem(lp, **configuration_parameters)
            else:
                status = local_solve_gurobi_problem(lp, **configuration_parameters)
        except:
            status = 'failed'

        if status in optimal_solution_strings:
            break
        # Prefer clearly optimal solutions but there may be other
        # acceptable results.  Back this up and come back to if we don't get
        # an optimal one.
        elif status in acceptable_solution_strings:
            acceptable_solution = the_solver.format_solution(lp, cobra_model)
            acceptable_solution_status = status 

    if status in optimal_solution_strings:
        the_solution = the_solver.format_solution(lp, cobra_model)
    elif type(acceptable_solution) != types.NoneType:
        the_solution = acceptable_solution
        the_solution.status = status
    elif status != 'infeasible':
        if error_reporting:
            print '%s failed: %s'%(solver, status)
        the_solution = the_solver.format_solution(lp, cobra_model)
        # Keep this in case local solvers have modified the status
        the_solution.status = status
    else:
        # Otherwise we have a solution that didn't fail
        # but isn't optimal or acceptable:
        # e.g. the solver didn't find a solution or
        # it is deemed not worthy, maybe due to quality issues
        # To be safe declare all these as infeasible
        # to force a new solution attempt
        the_solution = None
        the_solution.status = 'infeasible'
        if error_reporting:
            print '%s failed: infeasible'%(solver)          

    cobra_model.solution = the_solution
    return lp


def local_solve_cplex_problem(lp, **kwargs):
    """A performance tunable method for solving a problem
    

    """
    from cobra.solvers.parameters import parameter_mappings
    from cobra.solvers.cplex_solver import set_parameter, get_status

    parameter_mappings = parameter_mappings['cplex']

    # Update parameter settings if provided
    if kwargs:
        [set_parameter(lp, parameter_mappings[k], v)
         for k, v in kwargs.iteritems() if k in parameter_mappings]

    lp.solve()
    # Strictly enforce the bounds on the solution, if available (not available for LPs)
    if 'x_bound_error_max' in dir(lp.solution.get_quality_metrics()):
        if (lp.solution.get_quality_metrics()).x_bound_error_max <= kwargs['tolerance_feasibility']:
            status = get_status(lp)
        else:
            status = get_status(lp)
            if status not in acceptable_solution_strings:
                # Don't need to modify here
                status = status
            else:
                # Then we need to enforce x_bound_infeasible status
                # so this is not missed
                status = 'x_bound_infeasible'
                
    else:
        status = get_status(lp)
    return status


def local_solve_gurobi_problem(lp, **kwargs):
    """A performance tunable method for solving a problem
    

    """
    from cobra.solvers.parameters import parameter_mappings
    from cobra.solvers.gurobi_solver import set_parameter, get_status

    parameter_mappings = parameter_mappings['gurobi']

    # Update parameter settings if provided
    if kwargs:
        [set_parameter(lp, parameter_mappings[k], v)
         for k, v in kwargs.iteritems() if k in parameter_mappings]
    
    lp.optimize()
    status = get_status(lp)
    return status

def test_turnover_metabolite_constraints(cobra_model, **kwargs):
    """Test turnover metabolite reactions.  If verbose, inform
    on difficulties

    Arguments:
     cobra_model

    kwargs:
     solver
      solver, tolerance_optimality, tolerance_feasibility, tolerance_barrier    
     epsilon
     fraction_growth
     objective_sense
     the_FVA_reactions
     error_reporting
     number_of_processes
     FVA_verbose


    """
    # Set up
    solver, tolerance_optimality, tolerance_feasibility, tolerance_barrier, tolerance_integer = get_solver_parameters(**kwargs)

    if 'epsilon' in kwargs:
        epsilon = kwargs['epsilon'] 
    else:
        epsilon = tolerance_feasibility * 1.01

    if 'fraction_growth' in kwargs:
        fraction_growth = kwargs['fraction_growth'] 
    else:
        fraction_growth = 1.

    if 'objective_sense' in kwargs:
        objective_sense = kwargs['objective_sense'] 
    else:
        objective_sense = 'maximize'

    if 'the_FVA_reactions' in kwargs:
        the_FVA_reactions = kwargs['the_FVA_reactions'] 
    else:
        the_FVA_reactions = None

    if 'error_reporting' in kwargs:
        error_reporting = kwargs['error_reporting'] 
    else:
        error_reporting = None

    if 'number_of_processes' in kwargs:
        number_of_processes = kwargs['number_of_processes'] 
    else:
        number_of_processes = 1

    if 'FVA_verbose' in kwargs:
        FVA_verbose = kwargs['FVA_verbose'] 
    else:
        FVA_verbose = False
    
    if (len([x.id for x in cobra_model.reactions if x.id.startswith("IRRMILP_")]) > 0):
        FVA_verbose = True
                    
    FVA_for_turnover_metabolites = irreversible_flux_variability_analysis(cobra_model,
        fraction_of_optimum = fraction_growth,
        objective_sense = objective_sense,
        the_reactions = the_FVA_reactions,
        solver = solver,
        tolerance_optimality = tolerance_optimality,
        tolerance_feasibility = tolerance_feasibility,
        tolerance_barrier = tolerance_barrier,
        tolerance_integer = tolerance_integer,
        error_reporting = error_reporting,
        number_of_processes = number_of_processes,
        verbose = FVA_verbose)

    unsustainable_reactions = []
    # Flip back the lower bound on metabolites that pass this check.
    for reaction_id in FVA_for_turnover_metabolites:
        if FVA_for_turnover_metabolites[reaction_id]['maximum']:
            if (FVA_for_turnover_metabolites[reaction_id]['maximum'] < epsilon):
                unsustainable_reactions.append(reaction_id)     
            else:
                cobra_model.reactions.get_by_id(reaction_id).lower_bound = epsilon                              
        else:
            unsustainable_reactions.append(reaction_id)
    if len(unsustainable_reactions) > 0:
        print("The following metabolites cannot sustain the specified minimum flux, constraints removed:")
        for x in unsustainable_reactions:
            print(x)

import re
number_finder = re.compile("[\d]+\.?[\d]*")

class penalty_number:
    def __init__(self, value):
        self.str_value = str(value)
        self.value = float(value)
    
    def __add__(self, other):
        # addition is like OR
        return penalty_number(min(self.value, other.value))
    
    def __mul__(self, other):
        # multiplication is like AND
        return penalty_number(max(self.value, other.value))


def evaluate_penalty_string(penalty_string):
    """ penalty string will have:
        * 'or' statements which need to be converted to min
        * 'and' statements which need to be converted to max
    >>> evaluate_penalty("(1 and 2 and 3)")
    max(1, 2, 3)


    """
    # if there are no ands or ors, we are done
    penalty_string = penalty_string.lower()  # don't want to deal with cases
    
    if "and" not in penalty_string and "or" not in penalty_string:
        return eval(penalty_string)
    # we will replace AND and OR with addition and multiplication
    # equivalent to min/max
    penalty_string = penalty_string.replace("or", "+").replace("and", "*")
    # replace the numbers with the custom class which have overloaded AND/OR
    values = [penalty_number(i) for i in number_finder.findall(penalty_string)]
    values_strings = tuple("values[%i]" % i for i in range(len(values)))
    penalty_string = number_finder.sub("%s", penalty_string)
    penalty_string = penalty_string % values_strings
    return eval(penalty_string).value


def evaluate_penalties(cobra_model, new_cobra_model, expression_dict, expression_threshold):
    """ Evaluate penalties based on reactions or genes
    Arguments:
     cobra_model: original cobra model, needed to get original, unmodified reaction id's.
     new_cobra_model: cobra_model with final reaction id's; e.g. irreversible
     expression_dict: dictionary with gene or reaction keys and expression values
     expression_threshold: values below this will be penalized.

    Returns:
     penalties: a dictionary with reaction keys and penalty values
    
    """
    from cobra.core.Reaction import Reaction
    # Check whether we have reaction or gene penalties.  Default is gene.
    # Reaction penalties added to help in situations like flux minimization.
    # Rather than user specifying, nicer to have an automated check here.
    reaction_ids = [x.id for x in cobra_model.reactions]
    if ((sum([1 for x in expression_dict.keys() if x in reaction_ids]) >
         (0.5 * len(expression_dict.keys())))
        | (expression_threshold == None)):
        penalties = {}
        # Filter these against model reactions to be safe, this way we
        # can make sure to get the reverse reactions.  Assume 
        # that forward and reverse reactions are similarly penalized
        # if no reverse if given.
        for the_reaction in new_cobra_model.reactions:
            if the_reaction.id in expression_dict.keys():
                cur_penalty = expression_threshold - expression_dict[the_reaction.id]
                penalties[the_reaction.id] = max(0, cur_penalty)
                if 'reflection' in dir(the_reaction):
                    if type(the_reaction.reflection) == Reaction:
                        if the_reaction.reflection.id not in expression_dict.keys():
                            penalties[the_reaction.reflection.id] = max(0, cur_penalty)
            # Assign zero penalty to reactions not in the expression_dict
            elif the_reaction.id not in penalties.keys():
                penalties[the_reaction.id] = 0
    else:
        penalties = evaluate_gene_penalties(new_cobra_model, expression_dict, expression_threshold)
    return penalties


def evaluate_gene_penalties(new_cobra_model, expression_dict, threshold):
    """ Evaluate the penalty terms for the optimization problem.  The gene 
    expression data and acceptance threshold is used to determine whether 
    and how much penalty per unit flux to impose in the optimization for 
    reactions below the threshold.

    Arguments
     new_cobra_model:
     expression_dict:
     threshold:
    Returns:
     penalites: a dict with reaction id's as keys
      and penalty values

    
    """
    import re
    from copy import deepcopy
    # First we test for gene P/A calls.
    # Declare P/A based on cutoff
    gene_pa_dict = {}
    for cur_gene in expression_dict.keys():
        if cur_gene in new_cobra_model.genes:
            if expression_dict[cur_gene] > threshold:
                gene_pa_dict[cur_gene] = 1
            else:
                gene_pa_dict[cur_gene] = 0
    penalties = {}
    # We also keep a list to track if a reaction is inactivated
    # due to GE data
    reaction_inactivated = []
    for test_reaction in new_cobra_model.reactions:
        test_reaction = deepcopy(test_reaction)
        if len(test_reaction.gene_reaction_rule) == 0:
            penalties[test_reaction.id] = 0
            # We can move onto the next reaction if
            # there is no GRR
        else:
            # Otherwise, we will have to evaluate the GR
            the_gene_reaction_relation = deepcopy(test_reaction.gene_reaction_rule)
            
            # Maybe wasteful with regards to computational resources, 
            # but pool all the call info here for the moment
            reaction_gene_dict = {}
            for the_gene in test_reaction.get_gene():
                reaction_gene_dict[the_gene] = {}
                if the_gene.id in expression_dict:
                    reaction_gene_dict[the_gene]['measurement'] = expression_dict[the_gene.id]
                    reaction_gene_dict[the_gene]['present'] = gene_pa_dict[the_gene.id]
                else:
                    # If there is no expression data, we assume the gene is present, 
                    # an arguably conservative assumption
                    reaction_gene_dict[the_gene]['present'] = 1
                    # If the gene is present, this measurement is
                    # just a placeholder
                    reaction_gene_dict[the_gene]['measurement'] = threshold * 2
            # Now evaluate the reaction; if it's on, then we move on, but if we are
            # to prefer "inactivation," we will need to evaluate the penalty
            the_gene_reaction_relation = deepcopy(test_reaction.gene_reaction_rule)
            for the_gene in test_reaction.get_gene():
                the_gene_re = re.compile('(^|(?<=( |\()))%s(?=( |\)|$))'%re.escape(the_gene.id))
                if reaction_gene_dict[the_gene]['present'] == 0:
                    the_gene_reaction_relation = the_gene_re.sub('False', the_gene_reaction_relation)
                else:
                    the_gene_reaction_relation = the_gene_re.sub('True', the_gene_reaction_relation)
            if eval(the_gene_reaction_relation):
                # If protein expression evaluates as "true" set the penalty to zero
                penalties[test_reaction.id] = 0
            else:
                # Inactivation is the scenario requiring calculation.
                # We would like to determine the minimal total penalty associated
                # with turning the reaction on.  This can be done by
                # substituting & with a maximum operator in the GPR and applying a minimum operator 
                # in place of OR.  Gene present values are associated with zero 
                # penalty and gene false/absent values are associated with a nonzero 
                # penalty.
                the_gene_reaction_relation = deepcopy(test_reaction.gene_reaction_rule)
                for the_gene in test_reaction.get_gene():
                    the_gene_re = re.compile('(^|(?<=( |\()))%s(?=( |\)|$))'%re.escape(the_gene.id))
                    if reaction_gene_dict[the_gene]['present'] == 1:
                        the_gene_reaction_relation = the_gene_re.sub('0', the_gene_reaction_relation)
                    else:
                        the_gene_penalty = str(threshold - reaction_gene_dict[the_gene]['measurement'])
                        the_gene_reaction_relation = the_gene_re.sub(the_gene_penalty, the_gene_reaction_relation)
                penalties[test_reaction.id] = evaluate_penalty_string(the_gene_reaction_relation)
    return penalties


def check_solver(solver = 'cplex'):
    """ A simple function to verify the solver to employ.


    """
    from cobra import solvers
    
    available_solvers = solvers.solver_dict.keys()

    if solver not in available_solvers:
        if 'cplex' in available_solvers:
            print("Switched solver from " + solver + " to cplex.")            
            return 'cplex'
        # Gurobi should work well but it is not
        # tested as thoroughly as cplex
        elif 'gurobi' in available_solvers:
            print("Switched solver from " + solver + " to gurobi.")             
            return 'gurobi'
        # glpk works well for LPs but is usually not effective
        # for MILPs
        elif 'glpk' in available_solvers:
            print("Switched solver from " + solver + " to glpk.") 
            return 'glpk'
        # Otherwise, none available.
        # Also note Java_glpk not tested/supported
        else:
            print("No working solvers detected.")
            return None
    else:
        return solver


def get_solver_parameters(**kwargs):
    """ A simple function to return the solver tolerance parameters
    Note the tolerances will preferably be supplied, the defaults
    may not be consistent with those used previously

    """
    if 'solver' in kwargs:
        solver = check_solver(kwargs['solver'])
    else:
        solver = check_solver('glpk')

    if 'tolerance_optimality' not in kwargs:
        tolerance_optimality = 1E-7
    else:
        tolerance_optimality = kwargs['tolerance_optimality']    

    if 'tolerance_feasibility' not in kwargs:
        tolerance_feasibility = 1E-7
    else:
        tolerance_feasibility = kwargs['tolerance_feasibility']

    if 'tolerance_barrier' not in kwargs:
        tolerance_barrier = 1E-7 * 0.0001
    else:
        tolerance_barrier = kwargs['tolerance_barrier']

    if 'tolerance_integer' not in kwargs:
        tolerance_integer = integer_tolerances[solver]
    else:
        tolerance_integer = kwargs['tolerance_integer']
    
    return solver, tolerance_optimality, tolerance_feasibility, tolerance_barrier, tolerance_integer


def irreversible_reaction_knockout_analysis(cobra_model, **kwargs):
    """
    This function provides a method for testing for gene, reaction, or turnover metabolite
     requirements for an irreversible model 
    
    Model should have objective identified and constraints set before calling
     this function!

    Arguments:
     cobra_model: Preferably, an objective reaction should be identified.  If not, we will
      just pick one to get a solution for hot-starting.

    kwargs:
     the_reactions: a dict of {key: [reactions]}, where all reactions under a key are KO'd together.
                    default None will result in a test KO for model reactions

     next arguments are solver/toelrance parameters, described previously: solver, tolerance_optimality
                    tolerance_feasibility, tolerance_barrier, tolerance_integer.

     starting_model: a separate cobra_model is used for hot-starting if supplied.  This 
      can be useful if there are difficulties finding solutions for the model.  Note 
      the bounds for each reaction should lie within the bounds for the cobra_model.
      This might be useful to put starting guesses in a different area of the solution space.

    error_reporting = None,
    number_of_processes
    verbose

    Returns:
     variability_dict: a dictionary with the reaction id's as top-level keys,
      FVA ranges, and solver status.


    """
    # Need to copy the model because we're updating
    # reactions.
    from numpy import hstack, zeros, array
    from scipy import sparse
    from copy import deepcopy
    from cobra.core.Metabolite import Metabolite
    from cobra.core.Reaction import Reaction
    from time import time
    import types
    from math import floor

    continue_flag = True
    # INITIALIZATION

    if 'the_reactions' in kwargs:
        the_reactions = kwargs['the_reactions'] 
    else:
        the_reactions = None

    solver, tolerance_optimality, tolerance_feasibility, tolerance_barrier, tolerance_integer = get_solver_parameters(**kwargs)

    if 'starting_model' in kwargs:
        starting_model = kwargs['starting_model'] 
    else:
        starting_model = None

    if 'error_reporting' in kwargs:
        error_reporting = kwargs['error_reporting'] 
    else:
        error_reporting = None

    if 'number_of_processes' in kwargs:
        number_of_processes = kwargs['number_of_processes'] 
    else:
        number_of_processes = 1

    if 'verbose' in kwargs:
        verbose = kwargs['verbose'] 
    else:
        verbose = False

    if solver == 'cplex':
        from cobra.solvers.cplex_solver import get_status
    elif solver == 'gurobi':
        from cobra.solvers.gurobi_solver import get_status
    elif solver == 'glpk':
        from cobra.solvers.glpk_solver import get_status

    # Back up model so we can manipulate
    cobra_model = cobra_model.copy()

    # If no reactions were provided we will use the model reactions
    if not the_reactions:
        the_reactions = {x.id: [x] for x in cobra_model.reactions if not (x.id.startswith("TMS_")| x.id.startswith("EX_") | x.id.startswith("DM_") | x.id.startswith("IRRMILP_"))}
        for the_reaction_id in test_reactions.keys():
            the_reaction = cobra_model.reactions.get_by_id(the_reaction_id)
            if 'reflection' in dir(the_reaction):
                if type(the_reaction.reflection) == Reaction:
                    if ((not the_reaction.id.endswith("_reverse")) & (the_reaction.reflection.id.endswith("_reverse"))):
                        if (the_reaction.reflection.id in test_reactions.keys()) & (the_reaction.id in test_reactions.keys()):
                            test_reactions.pop(the_reaction.reflection.id)
    elif type(the_reactions) == dict:
        the_reactions_new = {}
        for the_key in the_reactions.keys():
            if hasattr(the_key, 'id'):
                the_new_key = the_key.id
            else:
                the_new_key = the_key
            the_list = the_reactions[the_key]
            if type(the_list) != list:
                the_list = [the_list]
            the_new_list = []
            for the_item in the_list:
                if hasattr(the_item, 'id'):
                    the_item = the_item.id
                the_new_list.append(the_item)
            the_new_list = map(cobra_model.reactions.get_by_id, the_new_list)
            the_reactions_new[the_new_key] = the_new_list
        the_reactions = the_reactions_new 
    elif type(the_reactions) == list:
        the_reactions_new = {}
        for the_key in the_reactions:
            if hasattr(the_key, 'id'):
                the_new_key = the_key.id
            else:
                the_new_key = the_key
            the_reactions_new[the_new_key] = [cobra_model.reactions.get_by_id(the_new_key)]
        the_reactions = the_reactions_new

    variability_dict = {}
    # Not yet supported
    if number_of_processes > 1:
        print("This FVA script developed for irreversible models is not set up for multiple processors.  Exiting...")
        continue_flag = False

    # Try to find a solution for hot starting.
    else:
        # If an objective is detected, convert it to
        # a constraint for FVA to use the fr of optimum
        nonzero_objective_list = [x.id for x in cobra_model.reactions if abs(x.objective_coefficient) > 0]
        nonzero_objective_coefficient = 1
        if (len(nonzero_objective_list) > 0):
            if (len(nonzero_objective_list) > 1):            
                # Numerically it is safer to leave the optimum
                # implicitly constrained rather than adding a constraint
                # that essentially repeats this.  This was observed to make
                # a difference for rare MILP optimizations that
                # could not be solved in cplex when the redundant constraint
                # was added in.    
                convert_objective_to_constraint(cobra_model,
                    objective_sense = 'maximize', 
                    fraction_of_optimum = 0,
                    copy_model = False,
                    new_reaction_name = "objective",
                    solver=solver,
                    tolerance_optimality = tolerance_optimality,
                    tolerance_feasibility = tolerance_feasibility,
                    tolerance_barrier = tolerance_barrier,
                    tolerance_integer = tolerance_integer,
                    bound_best_optimum = False)
                objective_reaction = cobra_model.reactions.get_by_id('objective')
                objective_reaction.objective_coefficient = nonzero_objective_coefficient
                nonzero_objective_list = [objective_reaction.id]
            else:
                # Don't need to add_sense_aware_bounds()
                objective_reaction = cobra_model.reactions.get_by_id(nonzero_objective_list[0])
                nonzero_objective_coefficient = objective_reaction.objective_coefficient
            best_solution = gim3e_optimize(cobra_model,
                solver=solver,
                objective_sense='maximize',
                tolerance_optimality=tolerance_optimality,
                tolerance_feasibility=tolerance_feasibility,
                tolerance_barrier = tolerance_barrier,
                tolerance_integer = tolerance_integer)            
        else:        
            # if so there's no objective
            # set to get a hot-start
            # solution for FVA runs           
            candidates = cobra_model.reactions.query("biomass")
            objective_reaction = None
            for test_reaction in candidates:
                if ((test_reaction.upper_bound != 0) | (test_reaction.lower_bound != 0)):
                    objective_reaction = test_reaction
            if objective_reaction == None:
                print('Unable to identify an objective for hot-starting FVA solutions.  Exiting.')
                continue_flag = False
                best_solution = None
            else:
                objective_reaction.objective_coefficient = 1.
                nonzero_objective_coefficient = 1
                nonzero_objective_list = [objective_reaction.id]
                best_solution = gim3e_optimize(cobra_model,
                    solver = solver,
                    objective_sense = 'maximize',
                    tolerance_optimality = tolerance_optimality,
                    tolerance_feasibility = tolerance_feasibility,
                    tolerance_barrier = tolerance_barrier,
                    tolerance_integer = tolerance_integer)

        # Some additional checks in case we want to start from an existing model
        start_from_main_model = True
        from cobra.core.Model import Model
        if type(starting_model) == Model:
            starting_model = deepcopy(starting_model)
            reactions_passed = 0
            reaction_ids = [x.id for x in cobra_model.reactions]
            starting_reaction_ids = [x.id for x in starting_model.reactions]
            if len(starting_reaction_ids) == len(reaction_ids):
                for the_reaction in reaction_ids:
                    if the_reaction in starting_reaction_ids:
                        if starting_model.reactions.get_by_id(the_reaction).lower_bound >= cobra_model.reactions.get_by_id(the_reaction).lower_bound:
                            upper_bound = cobra_model.reactions.get_by_id(the_reaction).upper_bound
                            starting_upper_bound = starting_model.reactions.get_by_id(the_reaction).upper_bound
                            if starting_upper_bound - upper_bound < 1:
                                if starting_upper_bound - upper_bound > 0:
                                    print("Warning, decreasing upper bound in starting_model for " + the_reaction + ".")
                                    starting_model.reactions.get_by_id(the_reaction).upper_bound = upper_bound
                                if starting_model.reactions.get_by_id(the_reaction).objective_coefficient == cobra_model.reactions.get_by_id(the_reaction).objective_coefficient:
                                    reactions_passed += 1
                                else:
                                    print("Mismatch objective " + the_reaction + ".")
                            else:
                                print("Mismatch upper bound " + the_reaction + ".")
                        else:
                            print("Mismatch lower bound " + the_reaction + ".")
            if reactions_passed == len(reaction_ids):
                metabolites_passed = 0
                metabolite_ids = [x.id for x in cobra_model.metabolites]
                starting_metabolite_ids = [x.id for x in starting_model.metabolites]
                if len(starting_metabolite_ids) == len(metabolite_ids):
                    for the_metabolite in metabolite_ids:
                        if the_reaction in starting_metabolite_ids:
                            metabolites_passed += 1
                if metabolites_passed == len(metabolite_ids):
                    start_from_main_model = False            

        # Identify the hot start
        if start_from_main_model:
            print("Ignoring any supplied alternate model and starting from main model.")
            starting_model = deepcopy(cobra_model)
        else:
            best_solution = gim3e_optimize(starting_model,
                    solver = solver,
                    objective_sense = 'maximize',
                    tolerance_optimality = tolerance_optimality,
                    tolerance_feasibility = tolerance_feasibility,
                    tolerance_barrier = tolerance_barrier,
                    tolerance_integer = tolerance_integer)
            if type(starting_model.solution) != types.NoneType:
                print("Hot starting from supplied alternate model.")
            else:
                starting_model == deepcopy(cobra_model)
                best_solution = gim3e_optimize(starting_model,
                    solver = solver,
                    objective_sense = 'maximize',
                    tolerance_optimality = tolerance_optimality,
                    tolerance_feasibility = tolerance_feasibility,
                    tolerance_barrier = tolerance_barrier,
                    tolerance_integer = tolerance_integer)


        # Back up the optimal solution to access later
        optimized_model = deepcopy(starting_model)
        
        # Finally, test the KO's
        if continue_flag:
            if verbose:
                start_time = time()
            the_sense_dict = {'maximize': 'maximum'}
            the_reaction_keys = the_reactions.keys()
            for the_reaction_key in the_reaction_keys:
                the_reaction_list = the_reactions[the_reaction_key]
                # This flag makes sure we try to find a solution with the best basis
                move_to_next = False
                try_to_find_solution = True
                # We try to start with the best solution
                # but will have to start froms scratch if we violate it
                # by altering the reverse reaction bounds.
                the_solution = best_solution                
                # Default to infeasible
                tmp_dict = {}
                for the_sense, the_description in the_sense_dict.iteritems():
                    tmp_dict[the_description] = None
                    tmp_dict[str(the_description)+"_"+"status"] = 'infeasible'
                current_reaction_list_upper_bound = []
                current_reaction_list_lower_bound = []
                current_reaction_list_reflection_upper_bound = []
                current_reaction_list_reflection_lower_bound = []
                move_to_next = False
                for the_reaction in the_reaction_list:
                    current_reaction_list_upper_bound.append(the_reaction.upper_bound)
                    current_reaction_list_lower_bound.append(the_reaction.lower_bound)
                    # If constraining flux to 0 would violate constraint, then we know this
                    # deactivation cannot result in a feasible solution.  Note the model
                    # is irreversible and all fluxes must be positive.
                    if the_reaction.lower_bound > 0:
                        the_solution = None
                        move_to_next = True
                        try_to_find_solution = False
                    the_reaction.upper_bound = 0
                    the_reaction.lower_bound = 0
                    if (optimized_model.solution.x_dict[the_reaction.id] != 0.):
                        # If we now violate a constraint we will need to start
                        # the solver from scratch.  Eg GLPK will pass
                        # off a solution in violation with an acceptable
                        # flag if we start in violation.
                        the_solution = None
                        # e.g. To explore variability in the forward reaction we constrain
                        # the reverse reaction to zero.
                    if 'reflection' in dir(the_reaction):
                        if type(the_reaction.reflection) == Reaction:
                            reflection_id = the_reaction.reflection.id
                            reflection = cobra_model.reactions.get_by_id(reflection_id)
                            current_reaction_list_reflection_upper_bound.append(cobra_model.reactions.get_by_id(reflection_id).upper_bound)
                            current_reaction_list_reflection_lower_bound.append(cobra_model.reactions.get_by_id(reflection_id).lower_bound)
                            if cobra_model.reactions.get_by_id(reflection_id).lower_bound > 0:
                                the_solution = None
                                move_to_next = True
                                try_to_find_solution = False
                            cobra_model.reactions.get_by_id(reflection_id).lower_bound = 0
                            cobra_model.reactions.get_by_id(reflection_id).upper_bound = 0
                            if (optimized_model.solution.x_dict[reflection_id] != 0.):
                                # If the old solution now violates a
                                # constraint we will need to start
                                # the solver from scratch.  Eg GLPK
                                # has an artifact such that GLPK will pass
                                # off a solution in violation with an acceptable
                                # flag if we start in violation.
                                the_solution = None
                        else:
                            current_reaction_list_reflection_upper_bound.append(0)
                            current_reaction_list_reflection_lower_bound.append(0)    
                    else:
                        current_reaction_list_reflection_upper_bound.append(0)
                        current_reaction_list_reflection_lower_bound.append(0)

                for the_sense, the_description in the_sense_dict.iteritems():
                    while ((not move_to_next) & (try_to_find_solution)):
                        gim3e_optimize(cobra_model,
                            solver = solver,
                            objective_sense = the_sense,
                            # Need to include new_objective
                            # here to update the_problem
                            # when hot-starting.
                            tolerance_optimality = tolerance_optimality,
                            tolerance_feasibility = tolerance_feasibility,
                            tolerance_barrier = tolerance_barrier,
                            tolerance_integer = tolerance_integer,
                            the_problem = the_solution,
                            error_reporting = error_reporting)                                
                        if type(cobra_model.solution) != types.NoneType:
                            tmp_dict[the_description] = cobra_model.solution.f
                            tmp_dict[str(the_description)+"_"+"status"] = cobra_model.solution.status
                            move_to_next = True
                        else:
                            # Give this one more try if we haven't reset the basis
                            if the_solution != None:
                                the_solution = None
                                move_to_next = False
                            else:
                                move_to_next = True
                                tmp_dict[the_description] = None
                                tmp_dict[str(the_description)+"_"+"status"] = 'failed'

                # Additional check for solutions?
                # Haven't implemented this one here as in irreversible_flux_variability_analysis
                # Since we only test for maximizing the objective with KO's

                variability_dict[the_reaction_key] = tmp_dict
                
                # Now reset the reaction bounds
                for the_index, the_reaction in enumerate(the_reaction_list):
                    #current_reaction_list_upper_bound.append(the_reaction.upper_bound)
                    #current_reaction_list_lower_bound.append(the_reaction.lower_bound)
                    the_reaction.upper_bound = current_reaction_list_upper_bound[the_index]
                    the_reaction.lower_bound = current_reaction_list_lower_bound[the_index]
                    if 'reflection' in dir(the_reaction):
                        if type(the_reaction.reflection) == Reaction:
                            (cobra_model.reactions.get_by_id(the_reaction.reflection.id).lower_bound) = current_reaction_list_reflection_lower_bound[the_index]
                            (cobra_model.reactions.get_by_id(the_reaction.reflection.id).upper_bound) = current_reaction_list_reflection_upper_bound[the_index]

                best_solution_status = get_status(best_solution)
                if (best_solution_status not in acceptable_solution_strings):
                    # Reset best_solution for future tries.  Note this list was converted to 1 element previously.
                    the_objective_reaction = nonzero_objective_list[0]
                    best_solution = gim3e_optimize(starting_model,
                        solver = solver,
                        objective_sense = 'maximize',
                        tolerance_optimality = tolerance_optimality,
                        tolerance_feasibility = tolerance_feasibility,
                        tolerance_barrier = tolerance_barrier,
                        tolerance_integer = tolerance_integer)
                            
                if verbose:
                    next_reaction_id = "...none.  All done."
                    the_index = the_reaction_keys.index(the_reaction_key)
                    pass_s = time() - start_time
                    remaining_s = pass_s * (len(the_reactions) - (the_index + 1)) / (the_index + 1)
                    remaining_h = floor((remaining_s)/3600)
                    remaining_m = floor(((remaining_s)/3600 - remaining_h) * 60)                    
                    if (the_index + 1) < len(the_reactions):
                        next_reaction_id = the_reaction_keys[the_index + 1]
                    print("Completed "+ str(the_index + 1 ) + " of " + str(len(the_reaction_keys)) + ", about to try " + next_reaction_id + ".  El: %0.0f s.  R: %0.0f hr %0.0f min." % (pass_s, remaining_h, remaining_m))                        

    return variability_dict


def convert_group_dict_to_reversible(the_reactions, cobra_model):
    """
    Take an FVA result dict, e.g. as returned by
    irreversible_reaction_knockout_analysis, and
    convert it to reversible.

    
    """
    from copy import deepcopy
    from cobra.core import Reaction
    
    reference_reactions = deepcopy(the_reactions)
    for the_key in reference_reactions.keys():
        for index, the_reaction in enumerate(reference_reactions[the_key]):
            if 'reflection' in dir(the_reaction):
                if type(the_reaction.reflection) == Reaction:
                    if the_reaction.id.endswith("_reverse"):                   
                        the_reactions[the_key].append(the_reaction.reflection)
                        pop_index = [x.id for x in the_reactions[the_key]].index(the_reaction.id)
                        the_reactions[the_key].pop(pop_index)
        if len(set([x.id for x in the_reactions[the_key]])) < len(the_reactions[the_key]):
            temp = list(set([x.id for x in the_reactions[the_key]]))
            new_reactions = []
            for x in the_reactions[the_key]:
                if x.id in temp:
                    new_reactions.append(x)
                    temp.pop(temp.index(x.id))
            the_reactions[the_key] = new_reactions
    reference_reactions = deepcopy(the_reactions)
    for the_key in reference_reactions.keys():
        if the_key.endswith("_reverse"):
            # This is a big flag, check and make sure the
            # forward, if it exists, isn't a copy...            
            if the_key in [x.id for x in cobra_model.reactions]:
                the_reaction = cobra_model.reactions.get_by_id(the_key)
                if 'reflection' in dir(the_reaction):
                    if type(the_reaction.reflection) == Reaction:
                        if not the_reaction.reflection.id.endswith("_reverse"):
                            if the_reaction.reflection.id in the_reactions.keys():
                                reverse_list = reference_reactions[the_key]
                                forward_list = reference_reactions[the_reaction.reflection.id]
                                if set([x.id for x in reverse_list]) == set([x.id for x in forward_list]):
                                    the_reactions.pop(the_reaction.id)

    return the_reactions
