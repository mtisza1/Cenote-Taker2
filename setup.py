#!/usr/bin/env python

from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="Cenote-Taker2",
    version="2.0.0",
    author="Michael J. Tisza",
    author_email="michael.tisza@gmail.com",
    description="Discovery and thorough annotation of viral genomes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mtisza1/Cenote-Taker2",
    packages=['Cenote-Taker2'],
    package_dir={'Cenote-Taker2': 'Cenote-Taker2'},
    include_package_data=True,
    license='Public Domain',
    zip_safe=False,

    install_requires=[
                      'python=3.5.5',
                      'prodigal=2.6.3',
                      'BWA=0.7.17',
                      'last=v1021',
                      'samtools=1.3',
                      'mummer=3.23',
                      'circlator=1.5.5',
                      'blast=2.9.0',
                      'hhsuite=3.2.0',
                      'bioawk=1.0',
                      'entrez-direct=13.3',
                      'krona=2.7.1',
                      'hmmer=3.3',
                      'bowtie2=2.3.5',
                      'trnascan-se=2.0.5',
                      'bbtools=37.62',
                      'tbl2asn=25.7',
                      'emboss=6.6.0'],

    entry_points={
        'console_scripts': [
            'Cenote-Taker2 = Cenote-Taker2.run_cenote-taker2:main'
        ]
    },

    classifiers=[
        'Programming Language :: Python :: 3 :: Only',
        "Operating System :: OS Independent",
        'Development Status :: 4 - Beta',
        'License :: Public Domain',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
        'Environment :: Console'
    ]
)