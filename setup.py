#!/usr/bin/python3
from setuptools import setup, find_packages

APP = ['abYdraw_compile_test.py']
DATA_FILES = []
OPTIONS = {'iconfile':'abydraw_icon.png.icns'}

setup(name='abYdraw',
    version='1.0',
    description="This is a programme designed to use our group's Antibody Markup Language (AbML) for describing bispecific antibody (BsAb) formats by either inputting an AbML descriptor string of a BsAb or by drawing a BsAb and outputting the its descriptor string. It is written in Python 3 and using standard packages TKinter to build the graphical interface in order to make it as accessible as possible.",
    author='James Sweet-Jones',
    author_email="ucbtjrs@ucl.ac.uk",
    app=APP,
    packages = find_packages(),
    options={'py2app': OPTIONS},
    setup_requires=['py2app'],
)