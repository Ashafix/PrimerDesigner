import os
import sys
import setuptools
from test import run_tests


class TestPrimerDesigner(setuptools.Command):
    description = "Automatically run the test suite for PrimerDesigner."
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        this_dir = os.getcwd()
        os.chdir("test")
        sys.path.insert(0, '')
        run_tests.main([])
        os.chdir(this_dir)


__version__ = "0.1"

with open("README.rst", "rb") as handle:
    readme_rst = handle.read().decode("ascii")


setuptools.setup(name='primer_designer',
                 version=__version__,
                 author='',
                 author_email='',
                 url='',
                 description='',
                 long_description=readme_rst,
                 download_url='',
                 classifiers=['Development Status :: 5 - Production/Stable',
                              'Intended Audience :: Developers',
                              'Intended Audience :: Science/Research',
                              'License :: Freely Distributable',
                              'Operating System :: OS Independent',
                              'Programming Language :: Python :: 3',
                              'Programming Language :: Python :: 3.4',
                              'Programming Language :: Python :: 3.5',
                              'Programming Language :: Python :: 3.6',
                              'Programming Language :: Python :: 3.6',
                              'Topic :: Scientific/Engineering',
                              'Topic :: Scientific/Engineering :: Bio-Informatics',
                              'Topic :: Software Development :: Libraries :: Python Modules',
                 ],
                 cmdclass={"test": TestPrimerDesigner},
                 )