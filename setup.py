from setuptools import setup, find_packages
from glob import glob

setup(
    name = "sequenceconcat",
    version = "0.0.1",
    scripts = glob( 'bin/*.py' ),
    package_dir = {'':'lib'},
    author = "Tyghe Vallard",
    author_email = "vallardt@gmail.com",
    description = "Handle concatting and splitting sequence files",
    license = "PSF",
    keywords = "",
    install_requires = ['biopython'],
    setup_requires = ['nose'],
    tests_require = ['nose','mock'],
)
