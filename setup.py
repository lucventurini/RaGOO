from setuptools import setup
import glob

scripts = glob.glob("*.py")

setup(
    name='RaGOO',
    version='1.1',
    description='A tool to order and orient genome assembly contigs via minimap2 alignments to a reference genome.',
    author='Michael Alonge',
    author_email='malonge11@gmail.com',
    packages=['ragoo_utilities'],
    package_dir={'ragoo_utilities': 'ragoo_utilities/'},
    install_requires=[
              'intervaltree',
              'numpy',
          ],
    scripts=scripts,
)