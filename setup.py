import ez_setup
ez_setup.use_setuptools()

from setuptools import setup, find_packages

exec(open('primavera/version.py').read()) # loads __version__

setup(
    name='primavera',
    version=__version__,
    author='Zulko',
    url='https://github.com/Edinburgh-Genome-Foundry/primavera',
    description='Primer selection + data analysis for DNA assembly validation',
    long_description=open('pypi-readme.rst').read(),
    license='see LICENSE.txt',
    keywords="DNA assembly sequencing primer design",
    packages=find_packages(exclude='docs'),
    install_requires=["numpy", "Biopython", "proglog", 'flametree', 'pandas',
                      'dna_features_viewer', 'dnachisel'])
