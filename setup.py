from setuptools import setup
from os.path import join, dirname

setup(
    name='OPBA2',
    version='1.0',
    packages=["orgpba2"],
    long_description=open(join(dirname(__file__), 'README.md')).read(),
)
