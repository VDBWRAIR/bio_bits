# Bootstrap setuptools
from ez_setup import use_setuptools
use_setuptools()

from setuptools import setup, find_packages
from glob import glob

setup(
    name = "bio_pieces",
    version = "0.0.7",
    packages = find_packages(),
    author = "Tyghe Vallard",
    author_email = "vallardt@gmail.com",
    description = "Handle concatting and splitting sequence files",
    license = "GPLv2",
    keywords = "biopython split fasta concat",
    install_requires = [
        'biopython',
        'argparse',
    ],
    setup_requires = [
        'nose'
    ],
    tests_require = [
        'nose',
        'mock',
        'tempdir',
    ],
    entry_points = {
        'console_scripts': [
            'sequence_concat = bio_pieces.sequence_concat:main',
            'sequence_files_concat = bio_pieces.sequence_files_concat:main',
            'sequence_split = bio_pieces.sequence_split:main',
            'cat_sequences = bio_pieces.cat_sequences:main',
            'phyml_seqrename = bio_pieces.phyml_seqrename:main',
            'raxmlrunner = bio_pieces.raxmlrunner:main',
            'phymlrunner = bio_pieces.phymlrunner:main',
        ],
    },
)
