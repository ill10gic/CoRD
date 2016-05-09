from setuptools import setup, find_packages

requires = [
    'appnope>=0.1.0',
    'backports.shutil-get-terminal-size>=1.0.0',
    'click>=6.6',
    'configparser>=3.3.0.post2',
    'cookies>=2.2.1',
    'Cython>=0.24',
    'decorator>=4.0.9',
    'funcsigs>=1.0.2',
    'gnureadline>=6.3.3',
    'h5py>=2.5.0',
    'hs-restclient>=1.2.2',
    'netCDF4>=1.2.3.1',
    'scipy==0.17.0',
    'nose>=1.3.7',
    'numpy==1.11.0',
    'oauthlib>=1.1.1',
    'pandas>=0.17.1',
    'pathlib2>=2.1.0',
    'pbr>=1.9.1',
    'pexpect>=4.0.1',
    'pickleshare>=0.7.2',
    'ptyprocess>=0.5.1',
    'python-dateutil>=2.5.0',
    'pytz>=2015.7',
    'requests>=2.10.0',
    'requests-oauthlib>=0.6.1',
    'requests-toolbelt>=0.6.0',
    'simplegeneric>=0.8.1',
    'six>=1.10.0',
    'traitlets>=4.2.1',
    'xlrd>=0.9.4'
]

tests_require = [
    'responses>=0.5.1',
]

setup(
    name='cord',
    description='Coupled RipCAS-DFLOW for Vegetation and Hydrology of Streams',
    version='0.1',
    author_email='maturner01@gmail.com',
    author='Matthew Turner',
    license='BSD3',
    classifiers=[
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 2.7',
        'Environment :: Console',
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering',
        'Intended Audience :: Science/Research'
    ],
    packages=find_packages(),
    include_package_data=True,
    install_requires=requires,
    tests_require=tests_require,
    entry_points='''
        [console_scripts]
        cord=cord.scripts.cord:cli
    '''
)
