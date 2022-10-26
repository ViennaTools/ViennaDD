#!/usr/bin/python3

from setuptools import setup, find_packages

myPackages=find_packages(
        where='.',
        include=['@PROJECT_ROOT_NAME@']
    )

# TODO: Read long_description from README.md

setup(name="@PROJECT_ROOT_NAME@",
    version="@PROJECT_ROOT_VERSION@",
    author="TU Wien - Institute for Microelectronics",
    packages=myPackages,
    package_dir= {'@PROJECT_ROOT_NAME@': '@PROJECT_ROOT_NAME@'},
    package_data = {'': ['*.so.*', '*.so', '*.py']},
    description="Example for C++ projects wrapped to python modules.",
)

