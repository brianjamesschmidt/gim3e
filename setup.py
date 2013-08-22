from setuptools import setup, find_packages
setup(
    name = "gim3epy",
    version = "1.0.2",
    packages = find_packages(),
)

import ez_setup
ez_setup.use_setuptools()
from setuptools import setup, find_packages
from pdb import set_trace
__version = '1.0.2'
setup(
    name = "gim3epy",
    version = __version,
    packages = find_packages(),
    #scripts = [''],
    setup_requires = [],
    install_requires = ['cobra>=0.2.0'],
    extras_require = {},

    package_data = {'': ['*.txt', '*.html','LICENSE','README','data/*','examples/*py']},

    author = "Brian J Schmidt",
    author_email = "brianjschmidt@gmail.com",
    description = "GIM3Epy is a package for the metabolic model-guided analysis of metabolomics and transcriptomics data",
    license = "GPL V3.0",
    keywords = "metabolism metabolomics transcriptomics genome modeling",
    url = "http://opencobra.sourceforge.net",
    test_suite = "cobra.test.suite",
    long_description = "This package contains the Python distribution for GIM3E (Gene Inactivation Moderated by Metabolism, Metabolomics, and Expression), GIM3Epy. GIM3E is an algorithm that enables the development of condition-specific models based on a genome-scale model of metabolism, an objective function, transcriptomics, and cellular metabolomics data. GIM3E establishes metabolite utilization requirements with metabolomics data, uses model-paired transcriptomics data to find experimentally supported solutions, and also provides calculations of the turnover (production / consumption) flux of metabolites.",

    )
    
