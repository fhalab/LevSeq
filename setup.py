from setuptools import setup, find_packages, Command
import os
import re

class CreateInitFile(Command):
    description = 'create __init__.py file in barcoding directory'
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        init_path = os.path.join('levseq', 'barcoding', '__init__.py')
        if not os.path.exists(init_path):
            with open(init_path, 'w') as f:
                f.write('# This file is automatically created during installation\n')
        print(f"Created {init_path}")

def read_version():
    path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'levseq/__init__.py')
    with open(path, 'r') as fh:
        return re.search(r'__version__\s?=\s?[\'"](.+)[\'"]', fh.read()).group(1)


def read_author():
    path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'levseq/__init__.py')
    with open(path, 'r') as fh:
        return re.search(r'__author__\s?=\s?[\'"](.+)[\'"]', fh.read()).group(1)


def read_email():
    path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'levseq/__init__.py')
    with open(path, 'r') as fh:
        return re.search(r'__author_email__\s?=\s?[\'"](.+)[\'"]', fh.read()).group(1)


def read_git():
    path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'levseq/__init__.py')
    with open(path, 'r') as fh:
        return re.search(r'__url__\s?=\s?[\'"](.+)[\'"]', fh.read()).group(1)


def readme():
    with open('README.md') as f:
        return f.read()


setup(name='levseq',
      version=read_version(),
      description='',
      long_description=readme(),
      long_description_content_type='text/markdown',
      author=read_author(),
      author_email=read_email(),
      cmdclass={
          'create_init': CreateInitFile,
      },
      url=read_git(),
      license='GPL3',
      project_urls={
          "Bug Tracker": read_git(),
          "Documentation": read_git(),
          "Source Code": read_git()},
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Natural Language :: English',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
      ],
      keywords=['Nanopore', 'ONT', 'evSeq'],
      packages=['levseq'],
      include_package_data=True,
      package_data={
          'levseq.barcoding': [
              'minion_barcodes.fasta',
              'demultiplex',
              'demultiplex-arm64'
              'demultiplex-x86',
          ],
      },
      entry_points={
          'console_scripts': [
              'levseq = levseq.cmd:main'
          ]
      },
      install_requires=['Bio',
                        'biopython',
                        'fsspec',
                        'h5py',
                        'holoviews',
                        'jupyterlab',
                        'mappy',
                        'matplotlib',
                        'ninetysix',
                        'numpy',
                        'pandas',
                        'pybedtools',
                        'pycoQC',
                        'pyfaidx',
                        'pyparsing',
                        'pysam',
                        'scipy',
                        'sciutil',
                        'seaborn',
                        'scikit-learn',
                        'statsmodels',
                        'tqdm'],
      python_requires='>=3.8',
      data_files=[("", ["LICENSE"])]
      )
