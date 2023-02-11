from setuptools import setup

setup(name='unTimely_Catalog_explorer',
      version='1.1.0',
      description='A search and visualization tool for the unTimely Catalog',
      url='https://github.com/fkiwy/unTimely_Catalog_explorer',
      author='Frank Kiwy',
      author_email='frank.kiwy@outlook.com',
      license='MIT',
      py_modules=['unTimely_Catalog_tools'],
      install_requires=['numpy', 'matplotlib', 'pillow', 'requests', 'astropy', 'reproject'])
