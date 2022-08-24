from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='dbg',
version='0.0.1',
description='Docking grid box generator for Autodock4 and Vina',
long_description=readme(),
long_description_type= 'text/markdown',
classifiers=[
    'Development Status :: 5 - Production/Stable',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3',
    'Operating System :: OS Independent'
],
url='',
author='UsmanJamshed',
author_email='jamshedu@berkeley.edu',
keywords='core package',
license='MIT',
packages=['dbg'],
install_requires=[],
include_package_data=True,
zip_safe=False)