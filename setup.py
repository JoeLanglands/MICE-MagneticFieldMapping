from setuptools import setup

setup(name='micemag',
      version='0.1',
      description='',
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
      
