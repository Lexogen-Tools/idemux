#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'importlib_resources; python_version < "3.7"',
    'dataclasses; python_version < "3.7"',
    'tqdm'
]

setup_requirements = ['pytest-runner',
                      'setuptools>=30.3.0',
                      'wheel',
                      'setuptools_scm']

test_requirements = ['pytest>=3',
                     "pytest-cov",
                     "pytest-pep8",
                     "coverage"]

setup(
    author="Falko Hofmann, Michael Moldaschl, Andreas Tuerk",
    author_email=('falko.hofmann@lexogen.com, michael.moldaschl@lexogen.com, '
                  'andreas.tuerk@lexogen.com'),
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8', 
        'Programming Language :: Python :: 3.9',
    ],
    description=("A Lexogen tool for demultiplexing and  index error correcting fastq "
                 "files. Works with Lexogen i7, i5 and i1 barcodes."),
    install_requires=requirements,
    long_description=readme + '\n\n' + history,
    long_description_content_type="text/x-rst",
    include_package_data=True,
    keywords='idemux',
    name='idemux',
    entry_points={
        'console_scripts': ['idemux=idemux.__main__:main'],
    },
    packages=find_packages(include=['idemux', 'idemux.*']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/lexogen-tools/idemux',
    zip_safe=False,
    use_scm_version={
        # duplicated config from pyproject.toml; keep in sync
        "write_to": "idemux/_version.py",
        "version_scheme": "post-release",
    },
    extras_require={
        'docs': [
            'sphinx'
        ],
        'dev': [
            'pytest',
            'pytest-cov',
            'coverage',
        ]
    }
)
