from setuptools import setup, find_packages

__version__ = "0.4"

setup(name='fragresp',
      author='Tobias Wulsdorf',
      author_email='tobias.wulsdorf@gmail.com',
      description='FragResp: A python program for automated decomposition of molecules \
and calculation of coherent resp charges.',
      version=__version__,
      license='MIT',
      platforms=['Linux'],
      zip_safe=False,
      packages=find_packages(),
      entry_points={
          'console_scripts':
              ['run_fragresp  = fragresp.scripts.run_fragresp:entry_point']},)