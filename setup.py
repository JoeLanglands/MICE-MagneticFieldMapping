from setuptools import setup

setup(name='micemag',
      version='0.1',
      description='Code for modelling the magnetic field of the superconducting magnets of the Muon Ionization Cooling Experiment',
      url='https://github.com/JoeLanglands/MICE-MagneticFieldMapping',
      author='Joe Langlands',
      author_email='j.langlands@sheffield.ac.uk',
      license='MIT',
      packages=['micemag'],
      install_requires=[
        'scipy',
        'numpy',
        'matplotlib',
        'iminiuit',
        ]
      )
      
