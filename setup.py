import ez_setup
ez_setup.use_setuptools()

from setuptools import setup, find_packages

exec(open('dna_sequencing_viewer/version.py').read()) # loads __version__

setup(name='dna_sequencing_viewer',
      version=__version__,
      author='Zulko',
    description='',
    long_description=open('README.rst').read(),
    license='see LICENSE.txt',
    keywords="",
    packages= find_packages(exclude='docs'))
