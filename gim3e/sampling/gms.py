# TODO
# This module is need of code cleanup/streamlining

from cobra import solvers
# Load these to help determine which solver solutions are OK
acceptable_solution_strings = ['optimal', 'MIP_optimal', 'optimal_tolerance']
# Removed 'x_bound_infeasible' as acceptable, this was causing problems when sampling.

# May want to include this as acceptable when running cplex
# if a small violation of the boundaries is OK: 'x_bound_infeasible'
optimal_solution_strings = ['optimal']
from  cobra import __version__ as cobra_version

from cobra.core.Object import Object
class sample_container(Object):
    def __init__(self, full_cobra_model, name=None):
        """An object for housing a cobra_model and additional fields for sampling

        
        """
        Object.__init__(self, name)
        del self.annotation
        del self.notes
        # TODO might want to clean up this class some more
        # but it works for now
        
        self.cobra_model_full = full_cobra_model

        # Will keep critical information on this solution object
        self.the_reaction_ids_full = [x.id for x in full_cobra_model.reactions]
        self.lb_full = [full_cobra_model.reactions.get_by_id(x).lower_bound for x in self.the_reaction_ids_full]
        self.ub_full = [full_cobra_model.reactions.get_by_id(x).upper_bound for x in self.the_reaction_ids_full]
        self.const_ind = [index for index, x in enumerate(self.the_reaction_ids_full) if self.ub_full[index] == self.lb_full[index]]
        self.const_values = [self.lb_full[x] for x in self.const_ind]
        self.reduced_model = ''
        self.the_reaction_ids_reduced = ''
        self.lb_reduced = ''
        self.ub_reduced = ''
        
    def make_reduced_model(self):
        # TODO
        # This is not tested/developed yet,
        # but ideally we reduce the dimensionality here
        # to remove inaccessible or fixed reactions
        from copy import deepcopy
        if type(self.reduced_model) == str:
            if cobra_version == '0.2.0':
                cobra_model = deepcopy(self.cobra_model_full)
            else:
                cobra_model = self.cobra_model_full.copy()
            for reaction_ind in self.const_ind:
                the_reaction_id = self.the_reaction_ids_full[reaction_ind]
                delete_me = cobra_model.reactions.get_by_id(the_reaction_id)
                delete_me.delete()
                cobra_model.reactions.remove(delete_me)
            self.cobra_model_reduced = cobra_model
            self.the_reaction_ids_reduced = [x.id for x in cobra_model.reactions]
            self.lb_reduced = [self.cobra_model_reduced.reactions.get_by_id(x).lower_bound for x in self.the_reaction_ids_reduced]
            self.ub_reduced = [self.cobra_model_reduced.reactions.get_by_id(x).upper_bound for x in self.the_reaction_ids_reduced]


def achr_sampler(sampling_object, **kwargs):
        """ Run the sampling.

        kwargs:
         'n_points': number of points to gather
         'point_stopping_condition': how long to let sampling run for each point, default is n_warmup_points.
                           Default is 2x the number of reversible axes.
                           FYI: we used 200 before in the COBRA toolbox.
                           If x < 1, this is treated as a mixing_frac to attempt to converge to
         'bias': not implemented
         'sample_reduced': not tested
         'solver_tolerance': some measure of the accuracy of the solver tolerance
                             used to generate the warmup points is needed 
         'max_time': execution will interrupt and points found will be returned.  
         'edge_buffer': how much buffer to leave around the boundaries.
                        default is just slightly larger than the solver tolerance .
         'n_cpus': maximum cpus to use, but parallel execution not yet implemented


        """
        from numpy import hstack, zeros, array, vstack, NAN, reshape
        from scipy import sparse, dot, r_, transpose
        from copy import deepcopy
        from cobra.core.Metabolite import Metabolite
        from cobra.core.Reaction import Reaction
        from time import time
        import types
        from math import floor, sqrt
        import random

        # First check all the kwargs
        # n_ponts: set this equal to the warmup points
        # bias: not implemented yet

        if 'sample_reduced' in kwargs:
            if kwargs['sample_reduced']:
                if cobra_version == '0.2.0':
                    cobra_model = deepcopy(sampling_object.cobra_model_reduced)                    
                else:
                    cobra_model = sampling_object.cobra_model_reduced.copy()
                sampling_ids = (sampling_object.the_reaction_ids_reduced)
                ub = sampling_object.ub_reduced
                lb = sampling_object.lb_reduced
            else:
                if cobra_version == '0.2.0':
                    cobra_model = deepcopy(sampling_object.cobra_model_full)
                else:
                    cobra_model = sampling_object.cobra_model_full.copy()
                sampling_ids = sampling_object.the_reaction_ids_full
                ub = sampling_object.ub_full
                lb = sampling_object.lb_full                
        else:
            if cobra_version == '0.2.0':
                cobra_model = deepcopy(sampling_object.cobra_model_full)
            else:
                cobra_model = sampling_object.cobra_model_full.copy()
            sampling_ids = sampling_object.the_reaction_ids_full
            ub = sampling_object.ub_full
            lb = sampling_object.lb_full

        if 'solver' in kwargs:
            solver = check_solver(kwargs['solver'])
        else:
            solver = check_solver('glpk')

        if 'solver_tolerance' not in kwargs:
            solver_tolerance = 1E-7
        else:
            solver_tolerance = kwargs['solver_tolerance']

        # Tolerance for reprojecting the current point into the null space
        null_space_tolerance = solver_tolerance * 1E-6
		
        if 'max_time' in kwargs:
            max_time = kwargs['max_time']
        else:
            max_time = 0

	# Steps must be larger than this to be effective
        if 'max_min_factor' in kwargs:
            max_min_factor = kwargs['max_min_factor']
        else:
            max_min_factor = 1E-6
        max_min_tol = solver_tolerance * max_min_factor
        
        # Ignore directions where the distance to the boundary
	# is small when computing allowed
        # steps.  Here, can be a bit more liberal.
        # Note the defaults have been altered from the
        # sampler for the MATLAB COBRA toolbox
        if 'tolerance_u_factor' in kwargs:
            tolerance_u_factor = kwargs['tolerance_u_factor']
        else:
            tolerance_u_factor = 1E-6
        tolerance_u = solver_tolerance * tolerance_u_factor

	# A check on the size of the 
        # components of the unit step vector
        # Small components are ignored when checking
        # for how far the point can be moved
        # Note the defaults have been altered from the
        # sampler for the MATLAB COBRA toolbox
        if 'tolerance_d_factor' in kwargs:
            tolerance_d_factor = kwargs['tolerance_d_factor']
        else:
            tolerance_d_factor = 1E-6
        tolerance_d = solver_tolerance * tolerance_d_factor

        # Final check buffer for the turnover metabolites
        edge_buffer_check_final = 0

        warmup_points = sampling_object.warmup_points
        # Number of warmup points
        n_rxns_in_model, n_warmup_points = warmup_points.shape


        # Compose the irreversible reaction matrix as a reversible 
        # format to better facilitate sampling calculations
        reaction_partner_dict = get_rev_optimization_partner_dict(cobra_model)

        if 'n_points' in kwargs:
            n_points = kwargs['n_points']
        else:
            n_points = round(len(reaction_partner_dict.keys()) * 2)

        if 'point_stopping_condition' in kwargs:
            point_stopping_condition = kwargs['point_stopping_condition']
            if point_stopping_condition >= 1:
                point_stopping_type = 'n_steps'
            else:
                point_stopping_condition = kwargs['point_stopping_condition']
                point_stopping_condition = max(abs(1 - point_stopping_condition), abs(point_stopping_condition))
                point_stopping_type = 'mix_frac'
                check_steps = round(len(reaction_partner_dict.keys()) * 0.5)
                stop_steps = 100 * check_steps
        else:
            point_stopping_condition = round(len(reaction_partner_dict.keys()) * 2)
            point_stopping_type = 'n_steps'
            # Otherwise, we could give each warmup point a shot.
            # Can in theory pare this down to ~200 as in previous ACHR
            # By default I tie this to the number of reactions.
            # In preliminary testing I had OK results with 50% of 
            # the # of warmup points.  1.5x is even safer.
            # May want to increase if mix_frac still looks bad

        print("Initializing matrices for sampling... ")        
        # Extract the reaction data from the warmup points
        irreversible_reaction_ids = []
        irreversible_reaction_indices = []
        tms_reaction_ids = []
        tms_reaction_indices = []
        irreversible_reaction_ub = []
        irreversible_reaction_lb = []
        tms_reaction_ub = []
        tms_reaction_lb = []        
        # penalty needs to be treated like the virtual metabolites
        penalty_reaction_id = 'penalty'
        for index, x in enumerate(sampling_ids):
            if ((not x.startswith("IRRMILP_")) & (not x.startswith("TMS_")) & (not x == penalty_reaction_id)):
                irreversible_reaction_ids.append(x)
                irreversible_reaction_indices.append(index)
                irreversible_reaction_ub.append(cobra_model.reactions.get_by_id(x).upper_bound)
                irreversible_reaction_lb.append(cobra_model.reactions.get_by_id(x).lower_bound)           
            elif not x.startswith("IRRMILP_"):
                tms_reaction_ids.append(x)
                tms_reaction_indices.append(index)
                tms_reaction_ub.append(cobra_model.reactions.get_by_id(x).upper_bound)
                tms_reaction_lb.append(cobra_model.reactions.get_by_id(x).lower_bound)
        tms_reaction_ub = array(tms_reaction_ub)
        tms_reaction_lb = array(tms_reaction_lb)

        forward_reaction_flag = zeros((len(irreversible_reaction_ids)))
        for index, the_reaction_id in enumerate(irreversible_reaction_ids):
            if not the_reaction_id.endswith('_reverse'):
                forward_reaction_flag[index] = 1

        # Make numpy arrays of the warmup data
        irreversible_reaction_warmup_points = warmup_points[irreversible_reaction_indices, :]
        tms_reaction_warmup_points = warmup_points[tms_reaction_indices, :]



        reversible_reaction_warmup_points, ub_rev, lb_rev, reversible_reaction_ids = convert_points_to_reversible(reaction_partner_dict, irreversible_reaction_ids, irreversible_reaction_ub, irreversible_reaction_lb, irreversible_reaction_warmup_points)

        convert_rev_to_irrev_array = create_rev_to_irrev_conversion_array(irreversible_reaction_ids, reversible_reaction_ids, reaction_partner_dict)
        calc_tms_from_irrev_array = create_tms_from_irrev_conversion_array(tms_reaction_ids, irreversible_reaction_ids, cobra_model)


        # Construct an S matrix for calculating the null space
        # Need a list of reactions and metabolites
        # Basically, need to covert the model to reversible
        # and remove TMS & penalty-type reactions
        reversible_model = convert_to_reversible(cobra_model)
        
        # Rows are metabolites and cols are reactions
        # these should match 1:1 to the model, just need to preserve the order
        reversible_model_reactions = [reversible_model.reactions.get_by_id(x) for x in reversible_reaction_ids]
        S_matrix = make_S_matrix(reversible_model_reactions, reversible_model.metabolites)
        N_rev = null(S_matrix)
                            
        ub_rev = array(ub_rev)
        lb_rev = array(lb_rev)
        center_rev_point = reversible_reaction_warmup_points.mean(1)
        previous_rev_point = deepcopy(center_rev_point)
        previous_center_rev_point = center_rev_point

        # Just throw this in here, will use it to build in dummy checks later
        continue_flag = True

        # Pre-allocate memory for the resulting points
        sampled_points = zeros((len(sampling_ids), n_points))
        sampled_points[:] = NAN
        initial_points = zeros((len(sampling_ids), n_points))
        initial_points[:] = NAN

        successful_steps = 0.
        from numpy import subtract as nsub
        from numpy import divide as ndiv
        from numpy import add as nadd
        from numpy import sum as nsum

        random.seed()
        if continue_flag:
            print("Starting to sample...")
            start_time = time()
            n_valid_points = 0
            while n_valid_points < n_points:
                attempted_steps_for_current_point = 0
                continue_moving_current_point = True
                # Create a random step size vector

                initial_rev_point = previous_rev_point + 0
                successful_steps_for_current_point = 0
                

                # attempted_steps_for_current_point < n_steps_per_point:
                while continue_moving_current_point == True: 
        
                    attempted_steps_for_current_point += 1 
                    # Pick a random warmup point
                    random_point_id = random.sample([x for x in range(0, n_warmup_points)], 1)[0]
                    random_rev_point = reversible_reaction_warmup_points[:, random_point_id]

                    # Get a direction from the center to the warmup
                    u = random_rev_point - center_rev_point
                    # Make this a unit length vector
                    u = u / sqrt(dot(u,u))

                    # Calculate starting distances to upper and lower bounds             
                    dist_ub = nsub(ub_rev, previous_rev_point)
                    dist_lb = nsub(previous_rev_point, lb_rev)

                    # Ignore if we start too close to a boundary
                    # this avoids issues with small directions
                    # Also ignore directions with no motion component
                    valid_dir_rev_ind = (dist_ub > tolerance_u) & (dist_lb > tolerance_u) & (u != 0)

                    # Figure out positive and negative directions
                    pos_direction = (u[valid_dir_rev_ind] > tolerance_d).nonzero()
                    neg_direction = (u[valid_dir_rev_ind] < -1 * tolerance_d).nonzero()

                    # Figure out all the possible maximum and minimum step sizes
                    ub_step_temp = ndiv(dist_ub[valid_dir_rev_ind], u[valid_dir_rev_ind])
                    lb_step_temp = ndiv(-1*dist_lb[valid_dir_rev_ind], u[valid_dir_rev_ind])
                    max_step_vec = r_[ub_step_temp[pos_direction], lb_step_temp[neg_direction]]
                    min_step_vec = r_[lb_step_temp[pos_direction], ub_step_temp[neg_direction]]

                    # Calculate the maximum steps in the direction of u until we hit a boundary
                    max_step = max_step_vec.min()
                    min_step = min_step_vec.max()

                    # Just move on: pick a new direction if the step is too small
		    # Otherwise, we will try to update
                    if not ((abs(min_step) < max_min_tol) & (abs(max_step) < max_min_tol)) | (min_step > max_step):
                        norm_step = random.random()
                        step_dist = (norm_step * (max_step - min_step)) + min_step
                        cur_rev_point = nadd(previous_rev_point, (step_dist * u))
                        # Reproject the point to avoid accumulating numerical errors
                        if attempted_steps_for_current_point % 25 == 0:
                            if (((abs(dot(S_matrix, cur_rev_point))).max()) > (null_space_tolerance)):
                                # Check only needed for magnitude, debugging, new should be smaller
                                cur_rev_point = dot(N_rev, dot(transpose(N_rev), cur_rev_point))

                        # Consider which points can be moved.  Just consider movable points 
			# to avoid complications from numerical errors in stable points
                        n_over = nsum(cur_rev_point[valid_dir_rev_ind] > ub_rev[valid_dir_rev_ind])
                        n_under = nsum(cur_rev_point[valid_dir_rev_ind] < lb_rev[valid_dir_rev_ind])
                        
                        if (n_over + n_under) < 1:
                             temp_rev_to_irrev = dot(convert_rev_to_irrev_array, cur_rev_point)
                             cur_irrev_point = (((temp_rev_to_irrev * forward_reaction_flag) > 0) * temp_rev_to_irrev) - (((temp_rev_to_irrev * (1 - forward_reaction_flag)) < 0) * temp_rev_to_irrev)
                             if len(tms_reaction_ids) > 0:
                                  cur_tms_point = dot(calc_tms_from_irrev_array, cur_irrev_point)
                                  tms_out_of_bounds = (cur_tms_point > (tms_reaction_ub + edge_buffer_check_final)) | (cur_tms_point < (tms_reaction_lb - edge_buffer_check_final))
                             else:
                                 tms_out_of_bounds = [0]

                             if nsum(tms_out_of_bounds) < 1:
                                 # recalculate the center point
                                 # print("Pass interim TMS check")
                                 previous_rev_point = cur_rev_point
                                 successful_steps += 1
                                 successful_steps_for_current_point += 1
                                 center_rev_point = ((n_warmup_points + successful_steps) * center_rev_point + cur_rev_point) / (n_warmup_points + successful_steps + 1)
                             # else:
                                 # print("Failed on interim TMS check")
                                 # Note a big reason for this may be the penalty if there is one.
                                 # So we will not update the point when this happens so we don't
                                 # wander into portions of the space that are not OK.
                
                        # else:
                             #print("Failed on interim border check")

                    if point_stopping_type == 'n_steps':
                        if successful_steps_for_current_point >= point_stopping_condition:
                            continue_moving_current_point = False
                    else:
                        # otherwise, point_stopping_type = mix_frac for current point
                        # to save time we won't include TMS here
                        if (successful_steps_for_current_point > 0) & (successful_steps_for_current_point % check_steps == 0):
                            if n_valid_points > 0:
                                # the_reversible_sampled[:, 0] = center_rev_point
                                the_reversible_sampled[:, n_valid_points + 1] = cur_rev_point
                                # the_reversible_initial[:, 0] = center_rev_point
                                the_reversible_initial[:, n_valid_points + 1] = initial_rev_point
                                
                                # the_reversible_sampled, the_reversible_sampled_list = convert_sampling_results_to_reversible(sampling_object, dont_keep_list = ["penalty"], keep_tms = False)
                                # the_reversible_initial, the_reversible_initial_list = convert_sampling_results_to_reversible(sampling_object, type = "initial", dont_keep_list = ["penalty"], keep_tms = False)
                                # the_reversible_sampled, the_reversible_sampled_list = concatenate_matrices_by_id(the_reversible_sampled[:, 0 : n_valid_points], the_reversible_sampled_list, reshape(previous_rev_point, (-1, 1)), reversible_model_reactions)
                                # the_reversible_initial, the_reversible_initial_list = concatenate_matrices_by_id(the_reversible_initial[:, 0 : n_valid_points], the_reversible_initial_list, reshape(initial_rev_point, (-1, 1)), reversible_model_reactions)

                            else:
                                # Row ids are in reversible_model_reaction_ids
                                the_reversible_sampled = zeros((len(reversible_model_reactions), (n_points + 1)))
                                the_reversible_sampled[:] = NAN
                                the_reversible_initial = zeros((len(reversible_model_reactions), (n_points + 1)))
                                the_reversible_initial[:] = NAN
                                the_reversible_sampled[:, 0] = center_rev_point
                                the_reversible_sampled[:, 1] = cur_rev_point
                                the_reversible_initial[:, 0] = center_rev_point
                                the_reversible_initial[:, 1] = initial_rev_point
    
                            # As a shortcut, check for points that are not movable in the current step
                            # rather than valid_dir_rev_ind, u == 0
                            cur_ind_fixed = [i for i, x in  enumerate(u == 0) if x == True]
                            # Do a quick check of the mixing fraction so far to make sure it looks OK
                            # const_ind_reversible = calculate_const_ind_reversible(sampling_object)
                            # Note the_reversible_sampled_list_temp == the_reversible_initial_list_temp
                            mix_frac = mix_fraction(the_reversible_sampled[:,0:(n_valid_points+2)], the_reversible_initial[:,0:(n_valid_points+2)], fixed = cur_ind_fixed)
                            # print(attempted_steps_for_current_point)
                            # print(mix_frac)
                            #    import pdb
                            #    pdb.set_trace()
                            if max((1 - mix_frac), mix_frac) <= point_stopping_condition:
                                continue_moving_current_point = False

                        elif successful_steps_for_current_point > stop_steps:
                            successful_steps_for_current_point = 0
                            attempted_steps_for_current_point = 0
                            cur_rev_point = initial_rev_point
                            print("Over step limit, reset.")
                                
                # Test the last successful point found during the steps
                cur_rev_point = previous_rev_point 
                # Add the resulting point to
                # the points found so far
                temp_rev_to_irrev = dot(convert_rev_to_irrev_array, cur_rev_point)
                cur_irrev_point = (((temp_rev_to_irrev * forward_reaction_flag) > 0) * temp_rev_to_irrev) - (((temp_rev_to_irrev * (1 - forward_reaction_flag)) < 0) * temp_rev_to_irrev)
                # Note that TMS reactions are effectively a non-convex constraint
                # for the reversible reaction solution space.
                # Better to run the sampling algorithm, screen for valid points,
                # the re-run.  In the future might want to
                # explore good algorithms for non-convex spaces.
                TMS_pass = True
                if len(tms_reaction_ids) > 0:
                    cur_tms_point = dot(calc_tms_from_irrev_array, cur_irrev_point)
                    # Check the TMS is still in bounds.
                    # Note that TMS requirements make the reversible
                    # axis space non-convex!
                    # This is an important check, this may fail since we didn't
                    # inform this when picking u.  May want to re-check the whole
                    # vector here
                    tms_out_of_bounds = (cur_tms_point > (tms_reaction_ub + edge_buffer_check_final)) | (cur_tms_point < (tms_reaction_lb - edge_buffer_check_final))
                    if nsum(tms_out_of_bounds) < 1:
                        previous_rev_point = cur_rev_point
                        # print("Pass final TMS check")
                    else:
                        # Check that one of the older points is not better
                        test_rev_point = previous_rev_point
                        TMS_pass = False
                        print("Failed on final TMS check")

                if TMS_pass:
                    # Initial values for the point, keep to calculate mixing
                    sampled_points[:, n_valid_points] = calculate_array_to_add(convert_rev_to_irrev_array, cur_rev_point, forward_reaction_flag, calc_tms_from_irrev_array, sampling_ids, irreversible_reaction_ids, tms_reaction_ids)
                    sampling_object.sampled_points = sampled_points
                    initial_points[:, n_valid_points] = calculate_array_to_add(convert_rev_to_irrev_array, initial_rev_point, forward_reaction_flag, calc_tms_from_irrev_array, sampling_ids, irreversible_reaction_ids, tms_reaction_ids)
                    sampling_object.initial_points = initial_points
                    
                    n_valid_points += 1
                    
                    # This will probably take a while to run so alert the user every so often
                    pass_s = time() - start_time
                    pass_h = floor((pass_s)/3600)
                    pass_m = floor(((pass_s)/3600 - pass_h) * 60)                         
                    # remaining_s = pass_s * (len(the_reactions) - (the_index + 1)) / (the_index + 1)
                    remaining_h = 0
                    remaining_m = 0
                    if max_time > 0:
                        remaining_s = (max_time - pass_s)
                        remaining_h = floor((remaining_s)/3600)
                        remaining_m = floor(((remaining_s)/3600 - remaining_h) * 60)
                        if pass_s > max_time:
                              return sampled_points
                        print("Valid point " + str(n_valid_points) + " of " + str(n_points) + ". El time is: %0.0f hr %0.0f min; max re is: %0.0f hr %0.0f min." % (pass_h, pass_m, remaining_h, remaining_m))
                    else:
                        print("Identified valid point " + str(n_valid_points) + " of " + str(n_points) + ". Elapsed time is: %0.0f hr %0.0f min." % (pass_h, pass_m))
        return sampled_points



def sampling_optimize(cobra_model, solver = 'cplex', error_reporting=True, **kwargs):
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
            # if tolerance is violated, report the result as 'infeasible'
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


def check_solver(solver = 'glpk'):
    """

    """
    
    available_solvers = solvers.solver_dict.keys()

    if solver not in available_solvers:
        if 'cplex' in available_solvers:
            print("Switched solver from " + solver + " to cplex.")            
            return 'cplex'
        # glpk works well for LPs but is usually not effective
        # for MILPs
        elif 'glpk' in available_solvers:
            print("Switched solver from " + solver + " to glpk.") 
            return 'glpk'
        # Gurobi should work well but it is not
        # tested as thoroughly as cplex
        elif 'gurobi' in available_solvers:
            print("Switched solver from " + solver + " to gurobi.")             
            return 'gurobi'
        # Otherwise, none available.
        # Also note Java_glpk not tested/supported
        else:
            print("No working solvers!  Exiting...")
            cobra_model.solution = None
            return None
    else:
        return solver


def get_rev_optimization_partner_dict(cobra_model):
    """ Return a simple dict to help with identifying the appropriate optimizations
    to perform on models that have been converted to an irreversible format

    Arguments:
     cobra_model

    Returns:
     reaction_partner_dict: each entry has keys:
      minimize: ID of reaction to minimize for the pair
      maximize: ID of reaction to maximize for the pair
      minimize_sense: optimization sense to effect minimization for the pair
      maximize_sense: optimization sense to effect maximization for the pair
      reverse: the reverse partner for the key
     

    """
    from cobra.core import Reaction

    the_reactions = [x.id for x in cobra_model.reactions if not x.id.startswith("IRRMILP_")]
    reaction_partner_dict = {}
    all_reactions_considered = []
    # Want to generate a list of rxn: forward, rev id's.  Then pick ones to generate warm up points.
    for the_reaction in the_reactions:
        if the_reaction not in all_reactions_considered:
            the_reaction_2 = ''
            the_reaction_1 = cobra_model.reactions.get_by_id(the_reaction)
            if 'reflection' in dir(the_reaction_1):
                if type(the_reaction_1.reflection) == Reaction:
                    the_reaction_2 = the_reaction_1.reflection.id
                    the_reaction_2 = cobra_model.reactions.get_by_id(the_reaction_2)
        if type(the_reaction_2) != str:
            if the_reaction_2.id.endswith('_reverse'):
                forward_reaction = the_reaction_1
                reverse_reaction = the_reaction_2
            else:
                forward_reaction = the_reaction_2
                reverse_reaction = the_reaction_1
            reaction_partner_dict[forward_reaction.id] = {}
            reaction_partner_dict[forward_reaction.id]['reverse'] = reverse_reaction.id
            if (reverse_reaction.lower_bound >= 0) & (forward_reaction.upper_bound == forward_reaction.lower_bound):
                # Here, we are constrained to operate with the reaction running in "reverse"
                # and can essentially ignore the forward reaction
                reaction_partner_dict[forward_reaction.id]['maximize'] = reverse_reaction.id
                reaction_partner_dict[forward_reaction.id]['maximize_sense'] = 'minimize'
                reaction_partner_dict[forward_reaction.id]['minimize'] = reverse_reaction.id
                reaction_partner_dict[forward_reaction.id]['minimize_sense'] = 'maximize'
            elif (forward_reaction.lower_bound >= 0) & (reverse_reaction.upper_bound == reverse_reaction.lower_bound):
                # Here, we are constrained to operate with the reaction running in "forward"
                # and can essentially ignore the reverse reaction
                reaction_partner_dict[forward_reaction.id]['maximize'] = forward_reaction.id
                reaction_partner_dict[forward_reaction.id]['maximize_sense'] = 'maximize'
                reaction_partner_dict[forward_reaction.id]['minimize'] = forward_reaction.id
                reaction_partner_dict[forward_reaction.id]['minimize_sense'] = 'minimize'
            else:
                # Otherwise we will assume both the forward and reverse reactions are
                # allowed
                if (reverse_reaction.lower_bound > 0) | (forward_reaction.lower_bound > 0):
                    print("Warning, potential numerical issues with the bounds for " + forward_reaction.id + " and " + reverse_reaction.id)
                reaction_partner_dict[forward_reaction.id]['maximize'] = forward_reaction.id
                reaction_partner_dict[forward_reaction.id]['maximize_sense'] = 'maximize'
                reaction_partner_dict[forward_reaction.id]['minimize'] = reverse_reaction.id
                reaction_partner_dict[forward_reaction.id]['minimize_sense'] = 'maximize'
        else:
            reaction_partner_dict[the_reaction_1.id] = {}
            reaction_partner_dict[the_reaction_1.id]['maximize'] = the_reaction_1.id
            reaction_partner_dict[the_reaction_1.id]['maximize_sense'] = 'maximize'                    
            reaction_partner_dict[the_reaction_1.id]['minimize'] = the_reaction_1.id
            reaction_partner_dict[the_reaction_1.id]['minimize_sense'] = 'minimize'
            # Dont assign this key if there is none: reaction_partner_dict[the_reaction_1.id]['reverse'] = ''
            all_reactions_considered.extend([the_reaction_1.id])          
    return reaction_partner_dict


def make_S_matrix(reactions, metabolites):
    """
    """
    # from scipy.sparse import dok_matrix
    from scipy import zeros
    coefficient_dictionary = {}
    S = zeros((len(metabolites), len(reactions)))
    for reaction_index, the_reaction in enumerate(reactions):
        for the_key, the_value in the_reaction._metabolites.items():
            metabolite_index = metabolites.index(the_key)
            S[metabolite_index, reaction_index] = the_value
    return S


def convert_to_reversible(cobra_model):
    """ Calculate a reversible model from an irreversible one

    
    """
    reactions_to_add = []
    from cobra.core.Reaction import Reaction
    from cobra.core import Metabolite
    from copy import deepcopy

    if 'penalty' in cobra_model.reactions:
        additional_reaction_id_deletion_list = ['penalty']
    else:
        additional_reaction_id_deletion_list = []
    if 'penalty' in cobra_model.metabolites:
        additional_metabolite_id_deletion_list = ['penalty']
    else:
        additional_metabolite_id_deletion_list = []
    deleted_list = []

    if cobra_version == '0.2.0':
        cobra_model = deepcopy(cobra_model)
    else:
        cobra_model = cobra_model.copy()      
    reaction_partner_dict = get_rev_optimization_partner_dict(cobra_model)
    non_irrmilp_reaction_ids = []
    for reaction_id in reaction_partner_dict.keys():
        non_irrmilp_reaction_ids.append(reaction_id)
        if 'reverse' in reaction_partner_dict[reaction_id].keys():
            non_irrmilp_reaction_ids.append(reaction_partner_dict[reaction_id]['reverse'])
    for the_metabolite_id in additional_metabolite_id_deletion_list:
        the_metabolite = cobra_model.metabolites.get_by_id(the_metabolite_id)
        the_metabolite.remove_from_model()
    for the_reaction_id in additional_reaction_id_deletion_list:
        the_reaction = cobra_model.reactions.get_by_id(the_reaction_id)
        the_reaction.delete()
        cobra_model.reactions.remove(the_reaction)
    model_reaction_ids = [x.id for x in cobra_model.reactions]
    for the_reaction_id in model_reaction_ids:
        the_reaction = cobra_model.reactions.get_by_id(the_reaction_id)
        if the_reaction.id not in non_irrmilp_reaction_ids:
            # Can get rid of these IRRMILP reactions
            the_reaction.delete()
            cobra_model.reactions.remove(the_reaction)
            deleted_list.append(the_reaction.id)
        else:
            if the_reaction.id.startswith('TMS_'):
                # Get rid of TMs as well
		# May need to add/take these out, had inconsistent results between models.
                for the_metabolite in the_reaction._metabolites.keys():
                    the_metabolite.remove_from_model()
                the_reaction.delete()
                cobra_model.reactions.remove(the_reaction)
                
    for the_reaction_id in reaction_partner_dict.keys():
        if 'reverse' in reaction_partner_dict[the_reaction_id].keys():
            forward_reaction = cobra_model.reactions.get_by_id(the_reaction_id)
            reverse_reaction = cobra_model.reactions.get_by_id(reaction_partner_dict[the_reaction_id]['reverse'])
            rev_ub = -1 * reverse_reaction.lower_bound
            rev_lb = -1 * reverse_reaction.upper_bound
            for_lb = forward_reaction.lower_bound
            for_ub = forward_reaction.upper_bound
            if rev_ub == rev_lb:
                # This would include the case where they are both 0
                cur_ub = for_ub
                cur_lb = for_lb
            elif for_ub == for_lb:
                cur_ub = rev_lb
                cur_lb = rev_ub
            elif for_ub > 0:
                cur_ub = for_ub
                if rev_lb < 0:
                    cur_lb = rev_lb
                else:
                    cur_lb = for_lb
            else:
                cur_ub = rev_ub
                cur_lb = rev_lb
            forward_reaction.upper_bound = cur_ub
            forward_reaction.lower_bound = cur_lb
            reverse_reaction.delete()
            cobra_model.reactions.remove(reverse_reaction)

    # clean up orphan metabolites
    the_metabolite_ids = [x.id for x in cobra_model.metabolites]
    for the_metabolite_id in the_metabolite_ids:
        the_metabolite = cobra_model.metabolites.get_by_id(the_metabolite_id)
        if len(the_metabolite._reaction) == 0:
            the_metabolite.remove_from_model()
        elif the_metabolite.id.startswith("IRRMILP_"):
            the_metabolite.remove_from_model()

    return cobra_model



def null(A, eps=1e-15):
    """ Quick function for calculating null space, see
    http://mail.scipy.org/pipermail/scipy-user/2005-June/004650.html

    """
    import scipy
    import scipy.linalg
    n, m = scipy.shape(A)
    if n > m :
        return scipy.transpose(null(scipy.transpose(A), eps))
        return null(scipy.transpose(A), eps)
    u, s, vh = scipy.linalg.svd(A)
    s=scipy.append(s,scipy.zeros(m))[0:m]
    null_mask = (s <= eps)
    null_space = scipy.compress(null_mask, vh, axis=0)
    return scipy.transpose(null_space)


def mix_fraction(sample1, sample2, **kwargs):
    """ Compares two sets of sampled points and determines how mixed
    they are.

    Arguments
     sample1, sample2   Ordered set of points, numpy arrays.  The points must be in
                       the same order otherwise it does not make sense.
    kwargs
     fixed (optional)   The directions which are fixed and are not expected (indices)

    Returns
     mix                the mix fraction.  Goes from 0 to 1 with 1 being
                        completely unmixed and .5 being essentially 
                        perfectly mixed.  


    """
    from numpy import min, isnan, median, ones, outer
    from scipy import ones
    
    if 'fixed' not in kwargs:
        fixed = []
    else:
        fixed = kwargs['fixed']

    # ignore NAN rows
    ignore_rows = isnan(min(sample1,1)) | isnan(min(sample2,1))
    if len(fixed) > 0:
        ignore_rows[fixed,:] = True
    keep_rows = ~ ignore_rows

    sample1_reduced = sample1[keep_rows,:]
    sample2_reduced = sample2[keep_rows,:]

    m1 = median(sample1_reduced, 1)
    LPproblem = median(sample2_reduced, 1)
    n_rxn_reduced, n_points = sample1_reduced.shape

    l1 = (sample1_reduced > (outer(m1, ones([1, n_points]))))
    eq1 = (sample1_reduced == outer(m1, ones([1, n_points])))
    l2 = (sample2_reduced > outer(LPproblem, ones([1, n_points])))
    eq2 = (sample2_reduced == outer(LPproblem, ones([1, n_points])))

    eqtotal = eq1 | eq2

    fr_mix = float(sum(sum((l1 == l2) & (~ eqtotal))))/float(l1.size-sum(sum(eqtotal)))

    return fr_mix

def extract_rev_matrix_from_irrev(sampling_structure, the_matrix):
    """ A function to extract a reversible matrix from an irreversible one.
    Useful for interpreting sampling results.

    """
    sampling_ids = sampling_object.the_reaction_ids_full
    cobra_model = sampling_object.cobra_model_full

    # Extract the reaction data from the warmup points
    irreversible_reaction_ids = []
    irreversible_reaction_indices = []
    tms_reaction_ids = []
    tms_reaction_indices = []
    irreversible_reaction_ub = []
    irreversible_reaction_lb = []
    tms_reaction_ub = []
    tms_reaction_lb = []
    
    # penalty needs to be treated like the virtual metabolites
    penalty_reaction_id = 'penalty'
    for index, x in enumerate(sampling_ids):
        if ((not x.startswith("IRRMILP_")) & (not x.startswith("TMS_")) & (not x == penalty_reaction_id)):
            irreversible_reaction_ids.append(x)
            irreversible_reaction_indices.append(index)
            irreversible_reaction_ub.append(cobra_model.reactions.get_by_id(x).upper_bound)
            irreversible_reaction_lb.append(cobra_model.reactions.get_by_id(x).lower_bound)           
        elif not x.startswith("IRRMILP_"):
            tms_reaction_ids.append(x)
            tms_reaction_indices.append(index)
            tms_reaction_ub.append(cobra_model.reactions.get_by_id(x).upper_bound)
            tms_reaction_lb.append(cobra_model.reactions.get_by_id(x).lower_bound)
    tms_reaction_ub = array(tms_reaction_ub)
    tms_reaction_lb = array(tms_reaction_lb)

    forward_reaction_flag = zeros((len(irreversible_reaction_ids)))
    for index, the_reaction_id in enumerate(irreversible_reaction_ids):
        if not the_reaction_id.endswith('_reverse'):
            forward_reaction_flag[index] = 1

    # Make numpy arrays from the warmup data
    irreversible_reaction_matrix = the_matrix[irreversible_reaction_indices, :]
    tms_reaction_matrix = the_matrix[tms_reaction_indices, :]
    tms_reaction_ids = [sampling_ids[i] for i in tms_reaction_indices]

    # Compose the irreversible reaction matrix as a reversible 
    # format to better facilitate sampling calculations
    reaction_partner_dict = get_rev_optimization_partner_dict(cobra_model)
    reversible_reaction_ids = [x for x in irreversible_reaction_ids if x in reaction_partner_dict.keys()]


    for row_index, the_reaction_id in enumerate(reversible_reaction_ids):
        index_1 = irreversible_reaction_ids.index(the_reaction_id)
        # Upper bound will be upper limit on maximize reaction.
        ub_f = irreversible_reaction_ub[index_1]
        lb_f = irreversible_reaction_lb[index_1]            
        if 'reverse' in reaction_partner_dict[the_reaction_id].keys():
            index_2 = irreversible_reaction_ids.index(reaction_partner_dict[the_reaction_id]['reverse'])
            net_vector = irreversible_reaction_warmup_points[index_1, :] - irreversible_reaction_warmup_points[index_2, :]
            lb_r = -1 * irreversible_reaction_ub[index_2]
            ub_r = -1 * irreversible_reaction_lb[index_2] 
            if ub_r == lb_r:
                # This would include the case where they are both 0
                cur_ub = ub_f
                cur_lb = lb_f
            elif ub_f == lb_f:
                cur_ub = lb_r
                cur_lb = ub_r
            elif ub_f > 0:
                cur_ub = ub_f
                if lb_r < 0:
                    cur_lb = lb_r
                else:
                    cur_lb = lb_f
            else:
                cur_ub = ub_r
                cur_lb = lb_r
            if row_index == 0:
                reversible_matrix_points = net_vector
                ub_rev = [cur_ub]
                lb_rev = [cur_lb]
            else:
                reversible_matrix_points = vstack([reversible_reaction_warmup_points, net_vector])
                ub_rev.append(cur_ub)
                lb_rev.append(cur_lb)                 
        else:
            if row_index == 0:
                reversible_matrix_points = irreversible_reaction_warmup_points[index_1][:]
                ub_rev = [ub_f]
                lb_rev = [lb_f]                    
            else:
                reversible_matrix_points = vstack([reversible_reaction_warmup_points, irreversible_reaction_warmup_points[index_1][:]])
                ub_rev.append(ub_f)
                lb_rev.append(lb_f)

    reversible_matrix_points = vstack([reversible_reaction_warmup_points, tms_reaction_matrix])
    reversible_reaction_ids += tms_reaction_ids
                
    return reversible_matrix_points, reversible_reaction_ids, ub_rev, ub_rev, tms_reaction_matrix, tms_reaction_ids, 


def save_sampling_object(the_sampling_object, filename, path = ""):
    """ Break apart and save sampling objects

    Arguments:
     the_sampling_object: object to save
     filename: effectively a prefix, don't include .pickle, .npy, etc..., these will be added
     path: directory to save to

    
    """
    import cPickle
    #import json
    from copy import deepcopy
    from numpy import save
    from numpy import load
    import os
    #the_sampling_object = deepcopy(the_sampling_object)
    if 'warmup_points' in dir(the_sampling_object):
        save((path + filename + "_warmup_point_matrix.npy"), the_sampling_object.warmup_points)
        # del the_sampling_object.warmup_points
    if 'sampled_points' in dir(the_sampling_object):
        save((path + filename + "_sampled_points.npy"), the_sampling_object.sampled_points)
        # del the_sampling_object.sampled_points
    if 'initial_points' in dir(the_sampling_object):
        save((path + filename + "_initial_points.npy"), the_sampling_object.initial_points)
        # del the_sampling_object.initial_points

    # Don't pickle the class object, this just creates difficulties
    # if you update the sampler
    # code.  Will update the class whenever loading.
    sampling_dict = {}

    if 'cobra_model_full' in dir(the_sampling_object):
        sampling_dict['cobra_model_full'] = the_sampling_object.cobra_model_full    
    if 'const_ind' in dir(the_sampling_object):
        sampling_dict['const_ind'] = the_sampling_object.const_ind    
    if 'const_values' in dir(the_sampling_object):
        sampling_dict['const_values'] = the_sampling_object.const_values       
    if 'lb_full' in dir(the_sampling_object):
        sampling_dict['lb_full'] = the_sampling_object.lb_full         
    if 'the_reaction_ids_full' in dir(the_sampling_object):
        sampling_dict['the_reaction_ids_full'] = the_sampling_object.the_reaction_ids_full
    if 'ub_full' in dir(the_sampling_object):
        sampling_dict['ub_full'] = the_sampling_object.ub_full
    if 'reduced_model' in dir(the_sampling_object):
        sampling_dict['reduced_model'] = the_sampling_object.reduced_model
    if 'the_reaction_ids_reduced' in dir(the_sampling_object):
        sampling_dict['the_reaction_ids_reduced'] = the_sampling_object.the_reaction_ids_reduced
    if 'lb_reduced' in dir(the_sampling_object):
        sampling_dict['lb_reduced'] = the_sampling_object.lb_reduced
    if 'ub_reduced' in dir(the_sampling_object):
        sampling_dict['ub_reduced'] = the_sampling_object.ub_reduced
        
    fp = open(path + filename + ".pickle", "wb")
    cPickle.dump(sampling_dict, fp)
    #json.dump(the_sampling_object, fp)
    fp.close()


def load_sampling_object(filename, path = ""):
    """ Load sampling objects from files

    Arguments:
     filename: effectively a prefix, don't include .pickle, .npy, etc..., these will be added
     dir: directory to load from

    
    """
    import cPickle
    import os
    from numpy import load
    
    fp = open(path + filename + ".pickle", "rb")
    the_sampling_dict = cPickle.load(fp)
    fp.close()

    the_sampling_object = sample_container(the_sampling_dict['cobra_model_full'])

    # Update in case there were 
    # changes made to the object

    if 'const_ind' in dir(the_sampling_dict):
        the_sampling_object.const_ind = sampling_dict['const_ind']
    if 'const_values' in dir(the_sampling_dict):
        the_sampling_object.const_values = sampling_dict['const_values']
    if 'lb_full' in dir(the_sampling_dict):
        the_sampling_object.lb_full = sampling_dict['lb_full']
    if 'the_reaction_ids_full' in dir(the_sampling_dict):
        the_sampling_object.the_reaction_ids_full = sampling_dict['the_reaction_ids_full']
    if 'ub_full' in dir(the_sampling_dict):
        the_sampling_object.ub_full = sampling_dict['ub_full']
    if 'reduced_model' in dir(the_sampling_dict):
        the_sampling_object.reduced_model = sampling_dict['reduced_model']
    if 'the_reaction_ids_reduced' in dir(the_sampling_dict):
        the_sampling_object.the_reaction_ids_reduced = sampling_dict['the_reaction_ids_reduced']
    if 'lb_reduced' in dir(the_sampling_dict):
        the_sampling_object.lb_reduced = sampling_dict['lb_reduced']
    if 'ub_reduced' in dir(the_sampling_dict):
        the_sampling_object.ub_reduced = sampling_dict['ub_reduced']
    
    if len(path) == 0:
        new_path = "."
    files = [f for f in os.listdir(new_path) if os.path.isfile(f)]
    if (filename + "_warmup_point_matrix.npy") in files:
        the_sampling_object.warmup_points = load(path + filename + "_warmup_point_matrix.npy")
    if (filename + "_sampled_points.npy") in files:
        the_sampling_object.sampled_points = load(path + filename + "_sampled_points.npy")
    if (filename + "_initial_points.npy") in files:
        the_sampling_object.initial_points = load(path + filename + "_initial_points.npy")        
    return the_sampling_object


def reduce_warmup_points(sampling_object, **kwargs):
    """ Check the warmup points to remove redundant vectors.
    
    kwargs:
    'solver_tolerance'
    'verbose'


    """
    import numpy
    from numpy import array, zeros, ones, nonzero
    from numpy import logical_and
    from copy import deepcopy

    if 'solver_tolerance' not in kwargs:
        solver_tolerance = 1E-7
    else:
        solver_tolerance = kwargs['solver_tolerance']

    if 'verbose' not in kwargs:
        verbose = True
    else:
        verbose = kwargs['verbose']        

    if 'sample_reduced' in kwargs:
        if kwargs['sample_reduced']:
            if cobra_version == '0.2.0':
                cobra_model = deepcopy(sampling_object.cobra_model_reduced)

            else:
                cobra_model = sampling_object.cobra_model_reduced.copy()
            sampling_ids = (sampling_object.the_reaction_ids_reduced)
            ub = sampling_object.ub_reduced
            lb = sampling_object.lb_reduced
        else:
            if cobra_version == '0.2.0':
                cobra_model = deepcopy(sampling_object.cobra_model_full)
            else:
                cobra_model = sampling_object.cobra_model_full.copy()
            sampling_ids = sampling_object.the_reaction_ids_full
            ub = sampling_object.ub_full
            lb = sampling_object.lb_full
    else:
        if cobra_version == '0.2.0':
            cobra_model = deepcopy(sampling_object.cobra_model_full)
        else:
            cobra_model = sampling_object.cobra_model_full.copy()
        sampling_ids = sampling_object.the_reaction_ids_full
        ub = sampling_object.ub_full
        lb = sampling_object.lb_full    

    warmup_points = sampling_object.warmup_points
    # Number of warmup points
    n_rxns_in_model, n_warmup_points = warmup_points.shape

    # Extract the reaction data from the warmup points
    irreversible_reaction_ids = []
    irreversible_reaction_indices = []
    tms_reaction_ids = []
    tms_reaction_indices = []
    irreversible_reaction_ub = []
    irreversible_reaction_lb = []
    tms_reaction_ub = []
    tms_reaction_lb = []        
    # penalty needs to be treated like the virtual metabolites
    penalty_reaction_id = 'penalty'
    for index, x in enumerate(sampling_ids):
        if ((not x.startswith("IRRMILP_")) & (not x.startswith("TMS_")) & (not x == penalty_reaction_id)):
            irreversible_reaction_ids.append(x)
            irreversible_reaction_indices.append(index)
            irreversible_reaction_ub.append(cobra_model.reactions.get_by_id(x).upper_bound)
            irreversible_reaction_lb.append(cobra_model.reactions.get_by_id(x).lower_bound)           
        elif not x.startswith("IRRMILP_"):
            tms_reaction_ids.append(x)
            tms_reaction_indices.append(index)
            tms_reaction_ub.append(cobra_model.reactions.get_by_id(x).upper_bound)
            tms_reaction_lb.append(cobra_model.reactions.get_by_id(x).lower_bound)
    tms_reaction_ub = array(tms_reaction_ub)
    tms_reaction_lb = array(tms_reaction_lb)

    forward_reaction_flag = zeros((len(irreversible_reaction_ids)))
    for index, the_reaction_id in enumerate(irreversible_reaction_ids):
        if not the_reaction_id.endswith('_reverse'):
            forward_reaction_flag[index] = 1

    # Make numpy arrays of the numerical components of the warmup data [no binary]
    test_matrix = sampling_object.warmup_points[irreversible_reaction_indices + tms_reaction_indices, :]
    delete_indices = []
    do_not_delete_indices = []

    from time import time
    ones_matrix = ones((1, n_warmup_points))
    from numpy import abs as nabs
    from numpy import subtract as nsubtract
    from numpy import sum as nsum
    from math import floor
    if verbose:
        start_time = time()
    for the_index in range((n_warmup_points - 1), -1, -1):
        if (the_index not in delete_indices) & (the_index not in do_not_delete_indices):
            the_difference_matrix = nabs(nsubtract((test_matrix[:, the_index, None] * ones_matrix), test_matrix))
            # If another point essentially repeats the index point then we can consider deleting the current point
            max_difference_per_vector = the_difference_matrix.max(0)
            if nsum(max_difference_per_vector > solver_tolerance) < (n_warmup_points - 1):
                # Indices that are bad, max difference is less than tolerance
                the_indices = nonzero(max_difference_per_vector < solver_tolerance)
                the_indices = the_indices[0].tolist()
                if len(set(do_not_delete_indices).intersection(set(the_indices))) == 0:
                    do_not_delete_indices.append(the_index)
                    # Remove the current index from the delete lis, get rid of the other replicates
                    the_indices.pop(the_indices.index(the_index))
                    delete_indices += the_indices
                else:
                    # Add in the indices not already destined for deletion
                    delete_indices += [x for x in the_indices if x not in do_not_delete_indices]
            if verbose:
                number_points_tested = n_warmup_points - the_index
                next_reaction_id = "...none.  All done."
                pass_s = time() - start_time
                remaining_s = (pass_s * (n_warmup_points) ) / number_points_tested
                remaining_h = floor((remaining_s)/3600)
                remaining_m = floor(((remaining_s)/3600 - remaining_h) * 60)
                if (the_index) > 0:
                    print("Completed "+ str(number_points_tested) + " of " + str(n_warmup_points) + " points.  El: %0.0f s.  R: %0.0f hr %0.0f min." % (pass_s, remaining_h, remaining_m))
    delete_indices = list(set(delete_indices))
    warmup_points = numpy.delete(warmup_points,delete_indices,1)
    return warmup_points


def convert_sampling_results_to_reversible(sampling_container, **kwargs):
    """ Take sampling results and convert them from a reversible to
    an irreversible format
	
	kwargs:
	 type: type of points to convert, "sampled," "initial," or "warmup"
         keep_tms: whether to keep TMS reactions.  default is true.
         dont_keep_list: list of additional reaction ids not to keep 

    Returns:
     the_converted_results
     converted_reaction_list

    """
    import types
    from numpy import zeros, NaN

    if "type" in kwargs:
		the_result_type = kwargs["type"]
    else:
		the_result_type = "sampled"

    if "dont_keep_list" in kwargs:
		dont_keep_list = kwargs["dont_keep_list"]
    else:
		dont_keep_list = []

    if "keep_tms" in kwargs:
		keep_tms = kwargs["keep_tms"]
    else:
		keep_tms = True
		
    # For now assume the full model will be converted to irreversible
    cobra_model = sampling_container.cobra_model_full
    sampled_point_names = sampling_container.the_reaction_ids_full

    if the_result_type == "warmup":
		sampled_matrix = sampling_container.warmup_points
    elif the_result_type == "initial":
		sampled_matrix = sampling_container.initial_points
    else:
		sampled_matrix = sampling_container.sampled_points
    
    the_raw_reactions_to_convert = [x for x in sampled_point_names if not x.startswith("IRRMILP_")]
    paired_reaction_dict = {}
    tested_reactions = []
    converted_reaction_list = []
    #reversible_partner_dict = get_rev_optimization_partner_dict(sampling_container.cobra_model_full)
    #reversible_reaction_ids = [x for x in sampling_container.the_reaction_ids_full if x in reversible_partner_dict.keys()]

    for the_reaction_id in the_raw_reactions_to_convert:
        the_reaction = cobra_model.reactions.get_by_id(the_reaction_id)
        if the_reaction not in tested_reactions:
            the_key = the_reaction.id
            the_values = [the_reaction_id]
        
            if 'reflection' in dir(the_reaction):
                if type(the_reaction.reflection) != types.NoneType:
                    if not the_reaction.reflection.id.endswith('_reverse'):
                        the_key = the_reaction.reflection.id
                    the_values.append(the_reaction.reflection.id)
            paired_reaction_dict[the_key] = the_values
            converted_reaction_list.append(the_key)
            tested_reactions += the_values

    n_original_reactions, n_points = sampled_matrix.shape
    the_converted_results = zeros((len(converted_reaction_list), n_points))
    the_converted_results[:] = NaN

    # Maintain the order of reactions set by the sampling_container
    # converted_reaction_list.sort()
    for the_final_index, the_key in enumerate(converted_reaction_list):
        the_forward_reaction = the_key
        the_forward_index = sampled_point_names.index(the_key)
        if len(paired_reaction_dict[the_forward_reaction]) == 2:
            the_reverse_reaction = paired_reaction_dict[the_key]
            the_reverse_reaction.pop(the_reverse_reaction.index(the_forward_reaction))
            the_reverse_reaction = the_reverse_reaction[0]
            the_reverse_index = sampled_point_names.index(the_reverse_reaction)
            net_flux = sampled_matrix[the_forward_index, :] - sampled_matrix[the_reverse_index, :]
        else:
            the_reverse_reaction = ''
            net_flux = sampled_matrix[the_forward_index, :]
        the_converted_results[the_final_index, :] = net_flux

    if keep_tms == False:
        discard_indices = []
        keep_indices = [x for x in range(0, len(converted_reaction_list))]
        for the_index in keep_indices:
            if converted_reaction_list[the_index].startswith("TMS_"):
                discard_indices.append(the_index)
        keep_indices = [i for i in keep_indices if i not in discard_indices]
        the_converted_results = the_converted_results[keep_indices, :]
        converted_reaction_list = [converted_reaction_list[i] for i in keep_indices]

    if len(dont_keep_list) > 0:
        discard_indices = [i for i, the_test_reaction in enumerate(converted_reaction_list) if the_test_reaction in (dont_keep_list)]
        keep_indices = [i for i in range(0, len(converted_reaction_list)) if i not in discard_indices]
        the_converted_results = the_converted_results[keep_indices, :]
        converted_reaction_list = [converted_reaction_list[i] for i in keep_indices]        

    return(the_converted_results, converted_reaction_list)
        
def update_bounds(sampling_object, **kwargs):
        """ Update the bounds for the sampling object based on a bounds_dict.
        If warmup_points are available, these will be checked to make sure
        they are included in the updated bounds.
        
        kwargs:
         'bounds_dict': a dictionary of {reaction_id: {'upper_bound': value, 'lower_bound': value}}
                        should not be reversible
         'update_reduced': not fully implemented

        """
        from numpy import hstack, zeros, array, concatenate
        from scipy import sparse
        from copy import deepcopy
        from cobra.core.Metabolite import Metabolite
        from cobra.core.Reaction import Reaction
        from time import time
        import types
        from math import floor
        from numpy import zeros
        import random
            
        if 'bounds_dict' in kwargs:
            bounds_dict = (kwargs['bounds_dict'])
        else:
            bounds_dict = {}

        enforce_warmup = True

        if 'warmup_points' in dir(sampling_object):
            use_warmup_points = True
            warmup_points = sampling_object.warmup_points
        else:
            use_warmup_points = False
            
        if 'update_reduced' in kwargs:
            if kwargs['update_reduced']:
                cobra_model = (sampling_object.cobra_model_reduced)                    
                sampling_ids = (sampling_object.the_reaction_ids_reduced)
                ub = sampling_object.ub_reduced
                lb = sampling_object.lb_reduced
            else:
                cobra_model = (sampling_object.cobra_model_full)
                sampling_ids = sampling_object.the_reaction_ids_full
                ub = sampling_object.ub_full
                lb = sampling_object.lb_full                
        else:
            cobra_model = deepcopy(sampling_object.cobra_model_full)
            sampling_ids = sampling_object.the_reaction_ids_full
            ub = sampling_object.ub_full
            lb = sampling_object.lb_full

        for the_reaction_id in bounds_dict:
            the_reaction_index = sampling_ids.index(the_reaction_id)
            if use_warmup_points:
                warmup_max = warmup_points[the_reaction_index,:].max()
                warmup_min = warmup_points[the_reaction_index,:].min()
                new_max = max(warmup_max, bounds_dict[the_reaction_id]['upper_bound'])
                new_min = min(warmup_min, bounds_dict[the_reaction_id]['lower_bound'])
            else:
                new_max = (bounds_dict[the_reaction_id]['upper_bound'])
                new_min = (bounds_dict[the_reaction_id]['lower_bound'])
            cobra_model.reactions.get_by_id(the_reaction_id).upper_bound = new_max
            cobra_model.reactions.get_by_id(the_reaction_id).lower_bound = new_min
            ub[the_reaction_index] = new_max
            lb[the_reaction_index] = new_min


def create_warmup_points(sampling_object, **kwargs):
        """ Create warm-up points for sampling.

        kwargs:
         'n_points': optional, will try to establish points to fully explore edges by FVA by default
          can try fewer or redundant samples if this is specified.  Generally don't want to specify
          this, let the sampler optimize the selection to make sure max/min for all reactions
          & TMs is included.
         'bias': not implemented yet
         'sample_reduced': not fully implemented yet.
         'solver'
         'solver_tolerance'
         'start_with_individual_reactions': boolean, default True, whether
         to completely randomize the warmup point optimization objectives or start with
         individual reactions.

        """
        from numpy import hstack, zeros, array, concatenate, ones, dot, NAN
        from scipy import sparse
        from copy import deepcopy
        from cobra.core.Metabolite import Metabolite
        from cobra.core.Reaction import Reaction
        from time import time
        import types
        from math import floor
        from numpy import random as nrandom
        import random

        if 'sample_reduced' in kwargs:
            if kwargs['sample_reduced']:
                if cobra_version == '0.2.0':
                    cobra_model = deepcopy(sampling_object.cobra_model_reduced)
                    
                else:
                    cobra_model = sampling_object.cobra_model_reduced.copy()
                sampling_ids = sampling_object.the_reaction_ids_reduced
            else:
                if cobra_version == '0.2.0':
                    cobra_model = deepcopy(sampling_object.cobra_model_full)
                    
                else:
                    cobra_model = sampling_object.cobra_model_full.copy()
                sampling_ids = sampling_object.the_reaction_ids_full
        else:
            if cobra_version == '0.2.0':
                cobra_model = deepcopy(sampling_object.cobra_model_full)
            else:
                cobra_model = sampling_object.cobra_model_full.copy()
            sampling_ids = sampling_object.the_reaction_ids_full
            
        if 'solver' in kwargs:
            solver = check_solver(kwargs['solver'])
        else:
            solver = check_solver('cplex')

        if 'solver_tolerance' not in kwargs:
            solver_tolerance = 1E-7
        else:
            solver_tolerance = kwargs['solver_tolerance']
        tolerance_optimality = solver_tolerance
        tolerance_feasibility = solver_tolerance
        tolerance_barrier = 0.0001 * solver_tolerance
        if solver == 'cplex':
            # We can set this to 0 in cplex
            tolerance_integer = 0
        else:
            tolerance_integer = 1E-9

        if 'start_with_individual_reactions' in kwargs:
            start_with_individual_reactions = (kwargs['start_with_individual_reactions'])
        else:
            start_with_individual_reactions = True


        the_bounds = None
        continue_flag = True
        print("Generating warm-up points...")
        
        # Don't include GIM3E MILP indicator type reactions
        reaction_partner_dict = get_rev_optimization_partner_dict(cobra_model)

        # now pick the start-up reactions.  first make a list of lists:
        # [rxn, min], [rxn, max]  then randomize the order &
        # we will march through in the random order
        test_list = []
        random.seed()
        for the_reaction in reaction_partner_dict.keys():
            maximize_reaction = reaction_partner_dict[the_reaction]['maximize']
            maximize_sense = reaction_partner_dict[the_reaction]['maximize_sense']
            if maximize_sense == 'maximize':
                maximize_coefficient = 1
            else:
                maximize_coefficient = -1
            minimize_reaction = reaction_partner_dict[the_reaction]['minimize']
            minimize_sense = reaction_partner_dict[the_reaction]['minimize_sense']
            if minimize_sense == 'maximize':
                minimize_coefficient = 1
            else:
                minimize_coefficient = -1

            if maximize_reaction != minimize_reaction:
                maximize_set = [[maximize_reaction, minimize_reaction], [maximize_coefficient, minimize_coefficient * -1]]
            else:
                maximize_set = [[maximize_reaction], [maximize_coefficient]]

            if minimize_reaction != maximize_reaction:
                minimize_set = [[minimize_reaction, maximize_reaction], [minimize_coefficient, maximize_coefficient * -1]]
            else:
                minimize_set = [[minimize_reaction], [minimize_coefficient]]

            # Couple the max/min, we will just keep a
            # Result if the while max-min pair works out
            test_list.append([maximize_set, minimize_set])

        the_sampling_guide = []
        if 'n_points' not in kwargs:
            print('Number of warmup points not specified, selecting one max/min pair for each of the axes.')
            n_points = len(test_list) * 2
            if start_with_individual_reactions:
                the_sampling_guide = test_list #random.sample(, n_points)
            
        else:
            n_points = kwargs['n_points']
            n_points = int(round(n_points/2.) * 2)
            if start_with_individual_reactions:
                the_sampling_guide = the_sampling_guide + random.sample(test_list, min(len(test_list), int(round(n_points/2.))))  
                          
        dim_x = len(sampling_ids)
        warmup_points = zeros((dim_x, 0))

        number_of_processes = 1
        variability_dict = {}
        # Not yet supported
        if number_of_processes > 1:
            print("This script developed for irreversible models is not set up for multiple processors.  Exiting...")
            continue_flag = False

        # Set up a hot start to speed things
        else:
            # If an objective is detected, will use it to generate a solution for hot-starting
            nonzero_objective_dict = {x.id: x.objective_coefficient for x in cobra_model.reactions if abs(x.objective_coefficient) > 0}         
            if (len(nonzero_objective_dict.keys()) > 0):
                # kwargs_to_pass = kwargs + {objective_sense: 'maximize'}
                best_solution = sampling_optimize(cobra_model,
                    solver = solver,
                    tolerance_optimality = tolerance_optimality,
                    tolerance_feasibility = tolerance_feasibility,
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
                    nonzero_objective_dict = {objective_reaction.id: 1}
                    best_solution = sampling_optimize(cobra_model,
                        solver = solver,
                        objective_sense = 'maximize',
                        new_objective = objective_reaction,
                        tolerance_optimality = tolerance_optimality,
                        tolerance_feasibility = tolerance_feasibility,
                        tolerance_barrier = tolerance_barrier,
                        tolerance_integer = tolerance_integer)

            # Back up the optimal solution to access later
            optimized_model = deepcopy(cobra_model)
            if cobra_version == '0.2.0':
                optimized_model = deepcopy(cobra_model)
            else:
                optimized_model = cobra_model.copy()   
            sampling_ids = sampling_object.the_reaction_ids_full
            for x in nonzero_objective_dict.keys():
                cobra_model.reactions.get_by_id(x).objective_coefficient = 0

        nrandom.seed()
        if continue_flag:
            verbose = True
            if verbose:
                start_time = time()

            test_reaction_id_list = reaction_partner_dict.keys()
            n_warmup_pts_found = 1
            number_trials_complete = 0
            # STEP THROUGH REACTIONS
            # Avoid duplicating warmup points
            while n_warmup_pts_found < n_points:
                # the_sample_reaction_id_list, the_sample_reactions_coefficients = the_sampling_guide[0]
                # the_reaction_coefficient_list = []

                warmup_point_pair = zeros((len(sampling_ids), 0))
                found_a_point_pair = False

                
                while not found_a_point_pair:
                    the_objective_dict = {}
                    for x in sampling_ids:
                        cobra_model.reactions.get_by_id(x).objective_coefficient = 0
                    # We try to start with the best solution
                    # but will have to start from scratch if we violate it
                    # by altering the reverse reaction bounds.
                    if len(the_sampling_guide) > 0:
                        for the_multiplier in [1, -1]:
                            the_objective_dict[the_multiplier] = {}
                            if the_multiplier == 1:
                                the_sample_reaction_id_list, the_sample_reactions_coefficients = the_sampling_guide[0][0]
                            else:
                                the_sample_reaction_id_list, the_sample_reactions_coefficients = the_sampling_guide[0][1]
                            for the_reaction_id in sampling_ids:
                                if the_reaction_id in the_sample_reaction_id_list:
                                    the_index = the_sample_reaction_id_list.index(the_reaction_id)
                                    the_objective_dict[the_multiplier][the_reaction_id] = the_sample_reactions_coefficients[the_index]
                                else:
                                    the_objective_dict[the_multiplier][the_reaction_id] = 0
                        the_sampling_guide.pop(0)
                    else:
                        warmup_coefficients = 2. * (nrandom.rand((len(test_reaction_id_list))) - 0.5)
                        for the_multiplier in [1, -1]:
                            the_objective_dict[the_multiplier] = {}
                            for the_reaction_index, the_reaction_id in enumerate(test_reaction_id_list):
                                the_coefficient = the_multiplier * warmup_coefficients[the_reaction_index]
                                if the_coefficient >= 0:
                                    the_maximize_reaction_id = reaction_partner_dict[the_reaction_id]['maximize']
                                    the_maximize_reaction_sense = reaction_partner_dict[the_reaction_id]['maximize_sense']
                                    the_reaction = cobra_model.reactions.get_by_id(the_maximize_reaction_id)
                                    if the_maximize_reaction_sense == 'maximize':
                                        the_objective_dict[the_multiplier][the_reaction.id] = the_coefficient
                                    if the_maximize_reaction_sense == 'minimize':
                                        the_objective_dict[the_multiplier][the_reaction.id] = -1. * the_coefficient
                                else:
                                    the_minimize_reaction_id = reaction_partner_dict[the_reaction_id]['minimize']
                                    the_minimize_reaction_sense = reaction_partner_dict[the_reaction_id]['minimize_sense']
                                    the_reaction = cobra_model.reactions.get_by_id(the_minimize_reaction_id)
                                    if the_minimize_reaction_sense == 'minimize':
                                        the_objective_dict[the_multiplier][the_reaction.id] = the_coefficient
                                    if the_minimize_reaction_sense == 'maximize':
                                        the_objective_dict[the_multiplier][the_reaction.id] = -1. * the_coefficient
                    for the_multiplier in [1, -1]:
                        first_try = True
                        found_a_point = False
                        the_solution = best_solution
                        
                        while not found_a_point:
                            for x in sampling_ids:
                                cobra_model.reactions.get_by_id(x).objective_coefficient = 0
                            for x in the_objective_dict[the_multiplier].keys():
                                cobra_model.reactions.get_by_id(x).objective_coefficient = the_objective_dict[the_multiplier][x]

                            the_result = sampling_optimize(cobra_model,
                                                       solver = solver,
                                                       objective_sense = 'maximize',
                                                       # Need to include new_objective
                                                       # here to update the_problem
                                                       # when hot-starting.
                                                       new_objective = the_objective_dict[the_multiplier],
                                                       the_problem = the_solution,
                                                       tolerance_optimality = tolerance_optimality,
                                                       tolerance_feasibility = tolerance_feasibility,
                                                       tolerance_barrier = tolerance_barrier,
                                                       tolerance_integer = tolerance_integer)

                            if type(cobra_model.solution.f) != types.NoneType:
                                if cobra_model.solution.status in optimal_solution_strings:
                                    found_a_point = True
                                    warmup_point_pair = concatenate((warmup_point_pair, zeros((len(sampling_ids),1))), axis=1)
                                    for storage_index, the_current_reaction_id in enumerate(sampling_ids):
                                         warmup_point_pair[storage_index, -1] = cobra_model.solution.x_dict[the_current_reaction_id]

                            if (type(cobra_model.solution.f) == types.NoneType):
                                # Reset the best_solution for future tries
                                for the_objective_reaction in sampling_ids:
                                    cobra_model.reactions.get_by_id(the_objective_reaction).objective_coefficient = 0
                                    if the_reaction_id in nonzero_objective_dict.keys():
                                        cobra_model.reactions.get_by_id(the_objective_reaction).objective_coefficient = nonzero_objective_dict[the_objective_reaction]
                                best_solution = sampling_optimize(cobra_model,
                                                              solver = solver,
                                                              the_problem = None,
                                                              new_objective = nonzero_objective_dict,
                                                              tolerance_optimality = tolerance_optimality,
                                                              tolerance_feasibility = tolerance_feasibility,
                                                              tolerance_barrier = tolerance_barrier,
                                    tolerance_integer = tolerance_integer)
                                for the_reaction_id in sampling_ids:
                                    the_reaction = cobra_model.reactions.get_by_id(the_reaction_id)
                                    the_reaction.objective_coefficient = 0
                                for x in the_objective_dict[the_multiplier].keys():
                                    cobra_model.reactions.get_by_id(x).objective_coefficient = the_objective_dict[the_multiplier][x]

                            if (not found_a_point) & (first_try == True):
                                # Try this one more time but force solving from scratch
                                the_solution = None
                                first_try = False
                            elif (not found_a_point) & (first_try == False):
                                # Give up and move on
                                found_a_point = True

                            if found_a_point & (the_multiplier == -1):
                                found_a_point_pair = True

                if warmup_point_pair.shape[1] == 2:
                    warmup_points = concatenate((warmup_points, warmup_point_pair), axis=1)
                    sampling_object.warmup_points = warmup_points
                    n_rows, n_warmup_pts_found = warmup_points.shape

                # Move onto the next

                if verbose:
                    pass_s = time() - start_time
                    if n_warmup_pts_found > 0:
                        remaining_s = (pass_s  / n_warmup_pts_found) * (n_points - n_warmup_pts_found)
                    else:
                        remaining_s = (pass_s * n_points) 
                    remaining_h = floor((remaining_s)/3600)
                    remaining_m = floor(((remaining_s)/3600 - remaining_h) * 60)
                    print("Found "+ str(n_warmup_pts_found) + " of " + str(n_points) + " warmup points.  El: %0.0f s.  R: %0.0f hr %0.0f min." % (pass_s, remaining_h, remaining_m))                  


            irreversible_reaction_ids = []
            irreversible_reaction_indices = []
            irreversible_reaction_ub = []
            irreversible_reaction_lb = []
            tms_reaction_ids = []
            tms_reaction_indices = []
            tms_reaction_ub = []
            tms_reaction_lb = []
            reaction_partner_dict = get_rev_optimization_partner_dict(cobra_model)
            
            penalty_reaction_id = 'penalty'
            for index, x in enumerate(sampling_ids):
                if ((not x.startswith("IRRMILP_")) & (not x.startswith("TMS_")) & (not x == penalty_reaction_id)):
                    irreversible_reaction_ids.append(x)
                    irreversible_reaction_indices.append(index)
                    irreversible_reaction_ub.append(cobra_model.reactions.get_by_id(x).upper_bound)
                    irreversible_reaction_lb.append(cobra_model.reactions.get_by_id(x).lower_bound)
                elif not x.startswith("IRRMILP_"):
                    tms_reaction_ids.append(x)
                    tms_reaction_indices.append(index)
                    tms_reaction_ub.append(cobra_model.reactions.get_by_id(x).upper_bound)
                    tms_reaction_lb.append(cobra_model.reactions.get_by_id(x).lower_bound)
            reversible_reaction_warmup_points, ub_rev, lb_rev, reversible_reaction_ids = convert_points_to_reversible(reaction_partner_dict, irreversible_reaction_ids, irreversible_reaction_ub, irreversible_reaction_lb, warmup_points)
            reversible_reaction_warmup_points_center = reversible_reaction_warmup_points.mean(1).reshape(len(reversible_reaction_ids),1)
            reversible_reaction_warmup_points = 0.33 * reversible_reaction_warmup_points + 0.67 * reversible_reaction_warmup_points_center * ones((1, n_points))
            convert_rev_to_irrev_array = create_rev_to_irrev_conversion_array(irreversible_reaction_ids, reversible_reaction_ids, reaction_partner_dict)
            calc_tms_from_irrev_array = create_tms_from_irrev_conversion_array(tms_reaction_ids, irreversible_reaction_ids, cobra_model)

            forward_reaction_flag = zeros((len(irreversible_reaction_ids)))
            for index, the_reaction_id in enumerate(irreversible_reaction_ids):
                if not the_reaction_id.endswith('_reverse'):
                    forward_reaction_flag[index] = 1

            warmup_points = zeros((len(sampling_ids), 0))

            for the_point_index in range(0, n_points):
                cur_rev_point = reversible_reaction_warmup_points[:, the_point_index]
                rev_to_irrev = dot(convert_rev_to_irrev_array, cur_rev_point)
                cur_irrev_point = (((rev_to_irrev * forward_reaction_flag) > 0) * rev_to_irrev) - (((rev_to_irrev * (1 - forward_reaction_flag)) < 0) * rev_to_irrev)
                cur_tms_point = dot(calc_tms_from_irrev_array, cur_irrev_point)
                array_to_add_warmup_points = zeros((len(sampling_ids),1))
                # Do it this way to set indicator rxns to nans since we aren't checking them
                for the_index, the_id in enumerate(sampling_ids):
                    the_value_sampled_points = NAN
                    the_value_initial_points = NAN
                    if the_id in irreversible_reaction_ids:
                        the_source_index = irreversible_reaction_ids.index(the_id)
                        the_value_warmup_points = cur_irrev_point[the_source_index]
                    elif the_id in tms_reaction_ids:
                        the_source_index = tms_reaction_ids.index(the_id)
                        the_value_warmup_points = cur_tms_point[the_source_index]
                array_to_add_warmup_points[the_index,0] = the_value_warmup_points
                warmup_points = concatenate((warmup_points, array_to_add_warmup_points), axis=1)



        return warmup_points



def convert_points_to_reversible(reaction_partner_dict, irreversible_reaction_ids, irreversible_reaction_ub, irreversible_reaction_lb, irreversible_reaction_warmup_points):
    """ Take a matrix of points from the irreversible model and convert to reversible format

    
    """
    from numpy import vstack
    reversible_reaction_ids = [x for x in irreversible_reaction_ids if x in reaction_partner_dict.keys()]
    for row_index, the_reaction_id in enumerate(reversible_reaction_ids):
        index_1 = irreversible_reaction_ids.index(the_reaction_id)
        # Upper bound will be upper limit on maximize reaction.
        ub_f = irreversible_reaction_ub[index_1]
        lb_f = irreversible_reaction_lb[index_1]            
        if 'reverse' in reaction_partner_dict[the_reaction_id].keys():
            index_2 = irreversible_reaction_ids.index(reaction_partner_dict[the_reaction_id]['reverse'])
            net_vector = irreversible_reaction_warmup_points[index_1, :] - irreversible_reaction_warmup_points[index_2, :]
            lb_r = -1 * irreversible_reaction_ub[index_2]
            ub_r = -1 * irreversible_reaction_lb[index_2] 
            if ub_r == lb_r:
                # This would include the case where they are both 0
                cur_ub = ub_f
                cur_lb = lb_f
            elif ub_f == lb_f:
                cur_ub = lb_r
                cur_lb = ub_r
            elif ub_f > 0:
                cur_ub = ub_f
                if lb_r < 0:
                    cur_lb = lb_r
                else:
                    cur_lb = lb_f
            else:
                cur_ub = ub_r
                cur_lb = lb_r
            if row_index == 0:
                reversible_reaction_warmup_points = net_vector
                ub_rev = [cur_ub]
                lb_rev = [cur_lb]
            else:
                reversible_reaction_warmup_points = vstack([reversible_reaction_warmup_points, net_vector])
                ub_rev.append(cur_ub)
                lb_rev.append(cur_lb)
        else:
            if row_index == 0:
                reversible_reaction_warmup_points = irreversible_reaction_warmup_points[index_1][:]
                ub_rev = [ub_f]
                lb_rev = [lb_f]                    
            else:
                reversible_reaction_warmup_points = vstack([reversible_reaction_warmup_points, irreversible_reaction_warmup_points[index_1][:]])
                ub_rev.append(ub_f)
                lb_rev.append(lb_f)
    return reversible_reaction_warmup_points, ub_rev, lb_rev, reversible_reaction_ids


def create_rev_to_irrev_conversion_array(irreversible_reaction_ids, reversible_reaction_ids, reaction_partner_dict):
    """ Array to multiple reversible array by to get irreversible back


    """
    from numpy import zeros
    # Create matrices that will be used to speed conversion from reversible
    # to irreversible format
    convert_rev_to_irrev_array = zeros((len(irreversible_reaction_ids), len(reversible_reaction_ids)))
    # Set the elements to ones
    # for the_irreversible_index in range(0, len(irreversible_reaction_ids)):
    for the_reversible_index in range(0, len(reversible_reaction_ids)):
        reversible_reaction_id = reversible_reaction_ids[the_reversible_index]
        # irreversible_reaction_id = irreversible_reaction_ids[the_irreversible_index]
        # First make sure these are paired
        the_dict = reaction_partner_dict[reversible_reaction_id]
        irr_index_1 = irreversible_reaction_ids.index(reversible_reaction_id)
        convert_rev_to_irrev_array[irr_index_1][the_reversible_index] = 1
        if 'reverse' in the_dict.keys():
            irr_index_2 = irreversible_reaction_ids.index(the_dict['reverse'])
            convert_rev_to_irrev_array[irr_index_2][the_reversible_index] = 1
    return convert_rev_to_irrev_array


def create_tms_from_irrev_conversion_array(tms_reaction_ids, irreversible_reaction_ids, cobra_model):
    """ We also want a matrix to convert the irreversible fluxes to
    TMS fluxes.


    """
    from numpy import zeros

    #  Repeat a similar trick, but don't need to
    # include the original irrev fluxes in the output vector
    calc_tms_from_irrev_array = zeros((len(tms_reaction_ids), len(irreversible_reaction_ids)))
    for the_tms_index in range(0, len(tms_reaction_ids)):
        the_tms_id = tms_reaction_ids[the_tms_index]
        the_tms = cobra_model.reactions.get_by_id(the_tms_id)
        # These should essentially be sink reactions with one element
        the_tms_coeff = -1. * [the_tms._metabolites[x] for x in the_tms._metabolites.keys()][0]
        the_vm_metabolite = [x for x in the_tms._metabolites.keys()][0]
        for the_reaction in the_vm_metabolite._reaction:
            if the_reaction.id != the_tms_id:
                irr_index = irreversible_reaction_ids.index(the_reaction.id)
                for the_metabolite in the_reaction._metabolites.keys():
                    if the_metabolite.id == the_vm_metabolite:
                        cur_vm_coeff = 1. * the_reaction._metabolites[the_metabolite]
                        calc_tms_from_irrev_array[the_tms_index][irr_index] = cur_vm_coeff / the_tms_coeff
    return calc_tms_from_irrev_array


def calculate_const_ind_reversible(sampling_object):
    """ Calcualte the constant indices for the reversible reactions.


    """
    the_reversible_points, the_reversible_list = convert_sampling_results_to_reversible(sampling_object, type = "warmup")
    reaction_partner_dict = get_rev_optimization_partner_dict(sampling_object.cobra_model_full)
    const_ind_reversible = sampling_object.const_ind
    the_irreversible_list = sampling_object.the_reaction_ids_full
    const_ind_reversible = []
    
    for the_reversible_index, the_reversible_id in enumerate(the_reversible_list):
        the_reversible_dict = reaction_partner_dict[the_reversible_id]
        if the_reversible_dict['maximize'] == the_reversible_dict['minimize']:
            the_irreversible_index = the_irreversible_list.index(the_reversible_id)
            if the_irreversible_index in const_ind_reversible:
                const_ind_reversible.append(the_irreversible_index)
        else:
            the_irreversible_index_1 = the_irreversible_list.index(the_reversible_dict['maximize'])
            the_irreversible_index_2 = the_irreversible_list.index(the_reversible_dict['minimize'])
            if the_irreversible_index_1 in const_ind_reversible:
                if the_irreversible_index_2 in const_ind_reversible:
                    const_ind_reversible.append(the_irreversible_index)

    return const_ind_reversible
                
        
        
def calculate_array_to_add(convert_rev_to_irrev_array, initial_rev_point, forward_reaction_flag, calc_tms_from_irrev_array, sampling_ids, irreversible_reaction_ids, tms_reaction_ids):
    """ This function pieces the results from sampling back together in a properly
    formatted new point for the array

    
    """

    from numpy import dot, zeros, NAN
    
    rev_to_irrev = dot(convert_rev_to_irrev_array, initial_rev_point)
    initial_irrev_point = (((rev_to_irrev * forward_reaction_flag) > 0) * rev_to_irrev) - (((rev_to_irrev * (1 - forward_reaction_flag)) < 0) * rev_to_irrev)
    cur_tms_initial_point = dot(calc_tms_from_irrev_array, initial_irrev_point)
    array_to_add_initial_points = zeros((len(sampling_ids)))
    # Do it this way to set indicator rxns to nans since we aren't checking them
    for the_index, the_id in enumerate(sampling_ids):
        the_value_sampled_points = NAN
        the_value_initial_points = NAN
        if the_id in irreversible_reaction_ids:
            the_source_index = irreversible_reaction_ids.index(the_id)
            the_value_initial_points = initial_irrev_point[the_source_index]
        elif the_id in tms_reaction_ids:
            the_source_index = tms_reaction_ids.index(the_id)
            the_value_initial_points = cur_tms_initial_point[the_source_index]
        array_to_add_initial_points[the_index] = the_value_initial_points


    return array_to_add_initial_points


def concatenate_matrices_by_id(the_matrix_1, the_matrix_1_row_id, the_matrix_2, the_matrix_2_row_id):
    """ This function joins two matrices based on lists of row ids
    
    """
    
    from numpy import reshape
    if the_matrix_1.ndim == 1:
        the_matrix_1 = reshape(the_matrix_1, (-1, 1))
    if the_matrix_2.ndim == 1:
        the_matrix_2 = reshape(the_matrix_2, (-1, 1))
    
    if the_matrix_1_row_id == the_matrix_2_row_id:
        from numpy import concatenate
        the_combined_list = the_matrix_1_row_id
        the_combined_matrix = concatenate((the_matrix_1, the_matrix_2), axis=1)
    else:
        from numpy import zeros, NAN, ix_
        the_combined_list = the_matrix_1_row_id
        the_combined_list += [x for x in the_matrix_2_row_id if x not in the_matrix_1_row_id]
        the_n_rows_1, the_n_cols_1 = the_matrix_1.shape
        the_n_rows_2, the_n_cols_2 = the_matrix_2.shape
        the_n_cols = the_n_cols_1 + the_n_cols_2
        the_n_rows = len(the_combined_list)
        the_combined_matrix = zeros((the_n_rows, the_n_cols))
        the_combined_matrix[:] = NAN
        the_indices_1 = [i for i, the_id in enumerate(the_combined_list) if the_id in the_matrix_1_row_id]
        the_combined_matrix[ix_(the_indices_1, [i for i in range(0, the_n_cols_1)])] = the_matrix_1
        the_indices_2 = [i for i, the_id in enumerate(the_combined_list) if the_id in the_matrix_2_row_id]
        the_combined_matrix[ix_(the_indices_2,[i for i in range(the_n_cols_1, the_n_cols)])] = the_matrix_2
    return the_combined_matrix, the_combined_list    
