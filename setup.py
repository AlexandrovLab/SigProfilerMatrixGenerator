from setuptools import setup, find_packages
import os
import shutil


#remove the dist folder first if exists
if os.path.exists("dist"):
	shutil.rmtree("dist")

def readme():
	with open('README.rst') as f:
		return(f.read())

VERSION = '1.0.14'

def write_version_py(filename='SigProfilerMatrixGenerator/version.py'):
	# Copied from numpy setup.py
	cnt = """
# THIS FILE IS GENERATED FROM SIGPROFILEMATRIXGENERATOR SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
	
	"""
	fh = open(filename, 'w')
	fh.write(cnt % {'version': VERSION,})
	fh.close()

write_version_py()

setup(name='SigProfilerMatrixGenerator',
		version=VERSION,
		description='SigProfiler matrix generator tool',
		url='',
		author='Erik Bergstrom',
		author_email='ebergstr@eng.ucsd.edu',
		license='UCSD',
		packages=find_packages(),#['SigProfilerMatrixGenerator'],
		install_requires =[
			"matplotlib>=2.2.2",
			"sigProfilerPlotting>=1.0.1",
			"statsmodels>=0.9.0",
			"scipy>=1.1.0",
			"pandas>=0.23.4",
			"numpy>=1.14.3"],
		include_package_data=True,
		zip_safe=False)
