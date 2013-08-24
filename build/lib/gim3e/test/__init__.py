from __future__ import with_statement, absolute_import
import sys
from os import name as __name
available_tests = ['unit_tests']

del __name

from os.path import abspath as __abspath
from os.path import join as __join
from os.path import split as __split
from os.path import sep as __sep
from cobra.manipulation import initialize_growth_medium

gim3e_directory = __abspath(__join(__split(__abspath(__file__))[0], ".."))
gim3e_location = __abspath(__join(gim3e_directory, ".."))
data_directory = gim3e_directory + "/data/"
gim3e_directory += '/core/'

salmonella_pickle = __join(data_directory, "salmonella_gem.pickle")
ecoli_sbml = __join(data_directory, "E_coli_core_M9.xml")
del __abspath, __join, __split, __sep

def create_test_model(test_pickle=salmonella_pickle):
    """Returns a cobra model for testing.  The default model is the
    version of the Salmonella enterica Typhimurium LT2 model published in
    Thiele et al. 2011 BMC Sys Bio 5:8, which has some updated metabolite
    KEGG id data for Schmidt et al. 2013 Bioinformatics

    test_pickle: The complete file name of a pickled cobra.Model or SBML XML
    file to be read.  We currently provide Salmonella enterica Typhimurium
    and Escherichia coli core models whose paths are stored in cobra.test.salmonella_pickle
    and cobra.test.ecoli_pickle, respectively.

    """
    from os import name as __name
    try:
        from cPickle import load
    except:
        from pickle import load

    with open(test_pickle, "rb") as infile:
        model = load(infile)
    initialize_growth_medium(model, 'LB')

    return model

def create_test_suite():
    """create a unittest.TestSuite with available tests"""
    from unittest import TestLoader, TestSuite
    loader = TestLoader()
    suite = TestSuite()
    for test_name in available_tests:
        exec("from . import " + test_name)
        suite.addTests(loader.loadTestsFromModule(eval(test_name)))
    return suite

suite = create_test_suite()

def test_all():
    """###running unit tests on gim3e###"""
    from unittest import TextTestRunner
    TextTestRunner(verbosity=2).run(create_test_suite())
