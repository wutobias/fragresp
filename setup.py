from setuptools import setup, find_packages

__version__ = "0.4"

setup(
    name='fragresp',
    author='Tobias Hüfner',
    author_email='tobias.huefner@biophys.mpg.de',
    description='FragResp: A python program for automated decomposition of molecules and calculation of coherent resp charges.',
    version=__version__,
    license='MIT',
    platforms=['Linux'],
    zip_safe=False,
    packages=find_packages(),
    entry_points={
        'console_scripts':
            [
                'run_fragresp=fragresp.scripts.run_fragresp:entry_point',
                'run_psi4=fragresp.scripts.run_psi4:entry_point'
            ]
        }
    )