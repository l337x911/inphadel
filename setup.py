from setuptools import setup, find_packages

setup(name='InPhaDel',
      version=1.0,
      author='Anand D. Patel',
      author_email='adp002@ucsd.edu',
      url='https://bitbucket.org/anand_patel/inphadel',
      packages=find_packages('src'),
      package_dir={'':'src'},
      package_data={'': ['src/models/*.pkl', '*.md',  
                              'src/abbrev_test/*.bed']},
      include_package_data=True,
      test_suite='svphase.test.predict.PredictionTest',
      zip_safe=True,
      install_requires=[
          # -*- Extra requirements: -*-
          "numpy>=1.6.0",
	  "pandas>=0.13.0",
	  "scikit-learn>=0.14.1"
      ],
      extras_require = {
        'MPL' : ["matplotlib>=1.1.1rc"]                 
      },
      entry_points = {
        'console_scripts': [ 
        'inphadel = svphase.inphadel:main',
        ]
      }
      )
