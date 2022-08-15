from setuptools import setup

CLASSIFIERS = [
    "Environment :: Console",
    "Environment :: MacOS X",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT license",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

setup(
 name='Strand Swaps',
 description="Find strand swaps of genes in a genbank file.",
 version='strand_swaps',
 author="Michael Roach",
 author_email="beardymcjohnface@gmail.com",
 py_modules=['strand_swaps'],
 install_requires=['Click>=7', 'snakemake>=6.10.0', 'pyyaml'],
 entry_points={
  'console_scripts': [
    'strand_swaps=strand_swaps.__main__:main'
  ]},
 include_package_data=True,
 package_data={"strand_swaps": ["Snakefile", "*.yml", "*.smk", "*.py", "*.R"]},
)
