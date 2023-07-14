from setuptools import setup, find_packages

setup(name='pubmedtools',
      version='0.0.0.2',
      author='Diogo de J. S. Machado',
      author_email='diogomachado.bioinfo@gmail.com',
      description=('Pubmed Tools (pubmedtools) package provides functions for'
                   'searching and retrieving articles from the PubMed database'
                   'using Biopython and NCBI Entrez Direct.'),
      packages=find_packages(),
      long_description_content_type='text/markdown',
      long_description=open('README.md').read(),
      zip_safe=False,
      install_requires=['pandas', 'biopython'],
      license = 'BSD-3-Clause',
      url='https://github.com/diogomachado-bioinfo/pubmedtools',
      )
