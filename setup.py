from setuptools import setup

with open('README.md') as f:
    long_description= f.read()

setup(name='dgbg',
version='0.0.1',
description='Docking grid box generator for Autodock4 and Vina',
long_description=long_description,
long_description_content_type= 'text/markdown',
classifiers=[
    'Development Status :: 5 - Production/Stable',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3',
    'Operating System :: OS Independent'
],
url='',
author='Usman Jamshed',
author_email='jamshedu@berkeley.edu',
keywords='core package',
license='MIT',
packages=['dgbg'],
install_requires=['pandas', 'numpy', 'mdtraj', 'nglview'],
include_package_data=True,
zip_safe=False,
download_url="https://github.com/ujamshed/dgbg/archive/refs/tags/0.0.1.tar.gz")