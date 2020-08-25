from setuptools import setup, find_packages

with open('README.md', 'r') as f:
    readme = f.read()

setup(
    name='tblastn_wrapper',
    version='0.1.0',
    description='A wrapper for tblastn that splits up queries from a file and runs tblastn on each individual query',
    long_description=readme,
    author='UNSW BINF6111 Team UCE',
    url='https://github.com/YasirKusay/tblastn_wrapper',
    packages=find_packages(exclude=('tests', 'docs'))
)
