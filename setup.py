from setuptools import setup, find_packages
import os

def readme():
	with open('README.rst') as f:
		return(f.read())


setup(name='SigProfilerMatrixGenerator',
		version='0.1.8',
		description='SigProfiler matrix generator tool',
		url='',
		author='Erik Bergstrom',
		author_email='ebergstr@eng.ucsd.edu',
		license='UCSD',
		packages=find_packages(),#['SigProfilerMatrixGenerator'],
		install_requires =[
			"matplotlib",
			"sigProfilerPlotting"],
		include_package_data=True,
		zip_safe=False)
