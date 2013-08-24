from __future__ import with_statement
import sys
from warnings import warn  # TODO - catch known warnings
from unittest import TestCase, TestLoader, TextTestRunner
from tempfile import gettempdir
from os import unlink
from os.path import join
from copy import deepcopy
from  cobra import __version__ as cobra_version
solver_tolerance = 1E-6
try:  #skipIf is not in python 2.6 / 2.5
    from unittest import skipIf
except:
    try:  # should we make unittest2 an absolute requirement and use skipIf
        from unittest2 import skipIf
    except:
        skipIf = None

# deal with absolute imports by adding the appropriate directory to the path
if __name__ == "__main__":
    sys.path.insert(0, "../..")
    from gim3e.test import data_directory, create_test_model
    from gim3e.test import ecoli_sbml as test_sbml_file
    from gim3e.test import salmonella_pickle as test_pickle
    from gim3e.core import gim3e
    sys.path.pop(0)
    assert 0
else:
    from . import data_directory, create_test_model, gim3e_directory
    from . import ecoli_sbml as test_sbml_file
    from . import salmonella_pickle as test_pickle
    sys.path.insert(0, gim3e_directory)
    import gim3e.core.gim3e as gim3e
    sys.path.pop(0)

# libraries which may or may not be installed
libraries = ["glpk", "gurobipy", "cplex"]
for library in libraries:
    try:
        exec("import %s" % library)
    except ImportError:
        exec("%s = None" % library)


class Gim3eTestCase(TestCase):
    def setUp(self):
        self.model = create_test_model(test_pickle)


class TestGim3eCore(Gim3eTestCase):
    """test core gim3e functions"""

    def test_solve_feasible(self):
        the_solver = gim3e.check_solver()
        tolerance_integer = gim3e.integer_tolerances[the_solver]
        if cobra_version == '0.2.0':
            test_model = deepcopy(self.model)
        else:
            test_model = self.model.copy()
        gim3e.gim3e_optimize(test_model, objective_sense = 'maximize', 
            the_problem=None, solver=the_solver,  
            error_reporting=None,
            tolerance_optimality = solver_tolerance, 
            tolerance_feasibility = solver_tolerance,
            tolerance_barrier=0.0001 * solver_tolerance,
            tolerance_integer = tolerance_integer)
        solution = test_model.solution        
        self.assertEqual(solution.status, "optimal")


    def test_solve_milp_feasible(self):
        the_solver = gim3e.check_solver()
        tolerance_integer = gim3e.integer_tolerances[the_solver]
        if cobra_version == '0.2.0':
            milp_model = deepcopy(self.model)
        else:
            milp_model = self.model.copy()
        gim3e.convert_to_irreversible_with_indicators(milp_model, mutually_exclusive_directionality_constraint = True)
        gim3e.gim3e_optimize(milp_model,objective_sense = 'maximize', 
            the_problem=None, solver=the_solver,  
            error_reporting=None,
            tolerance_optimality = solver_tolerance, 
            tolerance_feasibility = solver_tolerance,
            tolerance_barrier=0.0001 * solver_tolerance,
            tolerance_integer = tolerance_integer)
        solution = milp_model.solution        
        self.assertEqual(solution.status, "optimal")


    def test_optimize_turnover_metabolite(self):
        the_solver = gim3e.check_solver()
        tolerance_integer = gim3e.integer_tolerances[the_solver]
        if cobra_version == '0.2.0':
            milp_model = deepcopy(self.model)
        else:
            milp_model = self.model.copy()
        gim3e.convert_to_irreversible_with_indicators(milp_model, mutually_exclusive_directionality_constraint = False)
        epsilon = 1.01 * solver_tolerance
        gim3e.add_turnover_metabolites(milp_model, ['glc__D_c'], epsilon)
        for x in milp_model.reactions:
            x.objective_coefficient = 0
        milp_model.reactions.get_by_id('TMS_glc__D_c').objective_coefficient = 1
        gim3e.gim3e_optimize(milp_model, objective_sense = 'maximize', 
            the_problem=None, solver=the_solver,  
            error_reporting=None,
            tolerance_optimality = solver_tolerance, 
            tolerance_feasibility = solver_tolerance,
            tolerance_barrier=0.0001 * solver_tolerance,
            tolerance_integer = tolerance_integer)
        solution = milp_model.solution
        self.assertEqual(solution.status, "optimal")


    def test_convert_to_irreversible(self):
        from cobra.core.Reaction import Reaction
        reaction_id_to_check = [x.id for x in self.model.reactions]
        metabolite_id_to_check = [x.id for x in self.model.metabolites]
        gene_id_to_check = [x.id for x in self.model.genes]
        if cobra_version == '0.2.0':
            reversible_model = deepcopy(self.model)
        else:
            reversible_model = self.model.copy()
        # Since we won't solve for turnover eteabolites set mutually_exclusive_directionality_constraint = False
        gim3e.convert_to_irreversible_with_indicators(reversible_model, mutually_exclusive_directionality_constraint = False)
        for the_reaction_id in reaction_id_to_check:
            the_reaction = reversible_model.reactions.get_by_id(the_reaction_id)
            if 'reflection' in dir(the_reaction):
                if type(the_reaction.reflection) == Reaction:
                    the_reflection = the_reaction.reflection
                    # Check reactions point to the genes
                    gene_from_reaction = [x.id for x in the_reaction.get_gene()]
                    gene_from_reaction.sort()
                    gene_from_reflection = [x.id for x in the_reflection.get_gene()]
                    gene_from_reflection.sort()
                    self.assertEqual(gene_from_reaction, gene_from_reflection)
                    # Check that reactions point to the metabolites
                    metabolite_from_reaction = list(the_reaction.get_products() + the_reaction.get_reactants())
                    metabolite_from_reaction.sort()
                    metabolite_from_reflection = list(the_reflection.get_products() + the_reflection.get_reactants())
                    metabolite_from_reflection.sort()
                    self.assertEqual(metabolite_from_reaction, metabolite_from_reflection)
        for the_metabolite_id in metabolite_id_to_check:
            the_metabolite = reversible_model.metabolites.get_by_id(the_metabolite_id)
            the_reaction_id_list = [x.id for x in the_metabolite.get_reaction()]
            # Make sure forward and reverse reactions are both referenced
            for the_reaction_id in the_reaction_id_list:
                the_reaction = reversible_model.reactions.get_by_id(the_reaction_id)
                if 'reflection' in dir(the_reaction):
                    if type(the_reaction.reflection) == Reaction:
                        the_reflection = the_reaction.reflection
                        self.assertTrue(the_reflection.id in the_reaction_id_list)
        

    def test_convert_objective_to_constraint(self):
        the_solver = gim3e.check_solver()
        tolerance_integer = gim3e.integer_tolerances[the_solver]
        if cobra_version == '0.2.0':
            milp_model = deepcopy(self.model)
        else:
            milp_model = self.model.copy()
        gim3e.convert_to_irreversible_with_indicators(milp_model, mutually_exclusive_directionality_constraint = True)
        constraint_test_model = gim3e.convert_objective_to_constraint(milp_model,
                objective_sense = 'maximize', 
                fraction_of_optimum = 0.9,
                copy_model = True,
                bound_best_optimum = False,
                new_reaction_name = "test_objective",
                solver=the_solver,
                tolerance_optimality = solver_tolerance,
                tolerance_feasibility = solver_tolerance,
                tolerance_barrier = 0.0001 * solver_tolerance,
                tolerance_integer = tolerance_integer)
        constraint_test_model.reactions.get_by_id("test_objective").objective_coefficient = 1.
        gim3e.gim3e_optimize(milp_model, objective_sense = 'maximize', 
            the_problem=None, solver=the_solver,  
            error_reporting=None,
            tolerance_optimality = solver_tolerance, 
            tolerance_feasibility = solver_tolerance,
            tolerance_barrier=0.0001 * solver_tolerance,
            tolerance_integer = tolerance_integer)
        solution = milp_model.solution        
        self.assertEqual(solution.status, "optimal")


    def test_irreversible_flux_variability_analysis(self):
        the_solver = gim3e.check_solver()
        tolerance_integer = gim3e.integer_tolerances[the_solver]
        if cobra_version == '0.2.0':
            milp_model = deepcopy(self.model)
        else:
            milp_model = self.model.copy()
        gim3e.convert_to_irreversible_with_indicators(milp_model, mutually_exclusive_directionality_constraint = True)
        the_test_reaction_id = "biomass_iRR1083_metals"
        fva_dict = gim3e.irreversible_flux_variability_analysis(milp_model,
                fraction_of_optimum = 0.9,
                objective_sense= 'maximize',
                the_reactions = [the_test_reaction_id],
                solver = the_solver,
                tolerance_optimality = solver_tolerance,
                tolerance_feasibility = solver_tolerance,
                tolerance_barrier = 0.0001*solver_tolerance,
                tolerance_integer = tolerance_integer,
                error_reporting = None,
                number_of_processes = 1,
                verbose = False)
        self.assertEqual(fva_dict[the_test_reaction_id]['maximum_status'], "optimal")
        self.assertEqual(fva_dict[the_test_reaction_id]['minimum_status'], "optimal")


    def test_evaluate_penalties(self):
        the_solver = gim3e.check_solver()
        tolerance_integer = gim3e.integer_tolerances[the_solver]
        if cobra_version == '0.2.0':
            milp_model = deepcopy(self.model)
        else:
            milp_model = self.model.copy()
        gim3e.convert_to_irreversible_with_indicators(milp_model, mutually_exclusive_directionality_constraint = True)
        expression_dict = {x.id: 1. for x in self.model.genes}
        expression_dict.pop('s0001')
        penalties = gim3e.evaluate_penalties(self.model, milp_model, expression_dict, 2)
        for the_key in penalties.keys():
            if len(milp_model.reactions.get_by_id(the_key).gene_reaction_rule) > 0:
                if 's0001' not in (milp_model.reactions.get_by_id(the_key).gene_reaction_rule):
                    self.assertEqual(penalties[the_key], 1.)
                else:
                    self.assertEqual(penalties[the_key], 0.)
            else:
                self.assertEqual(penalties[the_key], 0.)


# make a test suite to run all of the tests
loader = TestLoader()
suite = loader.loadTestsFromModule(sys.modules[__name__])

def test_all():
    TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    test_all()
