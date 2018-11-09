import setuptools
from setuptools import setup, find_packages, Extension


from distutils.version import LooseVersion
if LooseVersion(setuptools.__version__) < LooseVersion('1.1'):
    print("Version detected:", LooseVersion(setuptools.__version__))
    raise ImportError(
        "the CGAT code collection requires setuptools 1.1 higher")


with open("README.md", "r") as fh:
    long_description = fh.read()

########################################################################
########################################################################
IS_OSX = sys.platform == 'darwin'

########################################################################
########################################################################
# collect CGAT version
sys.path.insert(0, "cgatshowcase")
import version

version = version.__version__

###############################################################
###############################################################
# Define dependencies
#
major, minor1, minor2, s, tmp = sys.version_info

if major < 3:
    raise SystemExit("""cgat-showcase requires Python 3 or later.""")

cgat_packages = find_packages()
cgat_package_dirs = {'cgatshowcase': 'cgatshowcase'}

##########################################################
##########################################################
# Classifiers
classifiers = """
Development Status :: 0 - Beta
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: MIT Licence
Programming Language :: Python :: 3
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""

setup(
    # package information
    name='cgatshowcase',
    version=version,
    description='cgatshowcase : the Computational Genomics Analysis Toolkit example pipeline/workflow',
    author='Adam Cribbs',
    author_email='adam.cribbs@imm.ox.ac.uk',
    license="MIT",
    platforms=["any"],
    keywords="computational genomics",
    long_description=long_description,
    classifiers=[_f for _f in classifiers.split("\n") if _f],
    url="https://github.com/cgat-developers/cgat-core",
    # package contents
    packages=cgat_packages,
    package_data={'cgatshowcase':['cgatshowcase/R/*.R']},
    include_package_dir=True,
    package_dir=cgat_package_dirs,
    include_package_data=True,
    # other options
    zip_safe=False,
    test_suite="tests",
)
