from setuptools import setup, find_packages
import os

def readme():
	with open('README.rst') as f:
		return(f.read())


setup(name='SigProfilerMatrixGenerator',
		version='0.1.21',
		description='SigProfiler matrix generator tool',
		url='',
		author='Erik Bergstrom',
		author_email='ebergstr@eng.ucsd.edu',
		license='UCSD',
		packages=find_packages(),#['SigProfilerMatrixGenerator'],
		install_requires =[
			"matplotlib==2.2.2",
			"sigProfilerPlotting==0.1.16",
			"statsmodels==0.9.0",
			"scipy==1.1.0",
			"pandas==0.23.4",
			"numpy==1.14.3"],
		include_package_data=True,
		zip_safe=False)
