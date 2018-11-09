# coding: utf-8
# Distributed under MIT License

import sys
import platform

from setuptools import setup, find_packages#, Extension
#from setuptools.command.build_ext import build_ext as _build_ext



extra_link_args = []
if sys.platform.startswith('win') and platform.machine().endswith('64'):
    extra_link_args.append('-Wl,--allow-multiple-definition')

long_desc = """
Hetero_Ct (Heterogeneous Catalysis) is a open-source Python library
for microkinetic modeling of heterogeneous catalyst reactions. 
"""

setup(
    name="hetero_ct",
    packages=find_packages(),
    version="0.0.1",
    #cmdclass={'build_ext': build_ext},
    setup_requires=['numpy>=1.14.3', 'setuptools>=18.0', 'cantera>=2.4.0'],
    install_requires=["numpy>=1.14.3", 'cantera>=2.4.0',
                      "monty>=0.9.6", "scipy>=1.0.1", 
                      "matplotlib>=1.5", "palettable>=2.1.1"],
    #extras_require={
    #    ':python_version == "2.7"': [
    #        'enum34',
    #    ],
    #    "provenance": ["pybtex"],
    #    "vis": ["vtk>=6.0.0"],
    #    },
    #package_data=None,
    author="Bharat Medasani",
    author_email="mbkumar@gmail.com",
    maintainer="Bharat Medasani",
    maintainer_email="mbkumar@gmail.com",
    #url=None,
    license="MIT",
    description="Hetero_Ct is a microkinetic modeling software targeting "
                "heterogeneous catalysis. It is based on Cantera.",
    long_description=long_desc,
    #keywords=["microkinetic modeling", "heterogeneous catalysis", "catalysis", "science",
    #          "project",  "surface"],
    classifiers=[
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Development Status :: 0 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Scientific/Engineering :: Chemical Engineering",
        "Topic :: Software Development :: Libraries :: Python Modules"
    ]#,
    #ext_modules=None,
    #entry_points=None
)
