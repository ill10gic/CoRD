from setuptools import setup, find_packages

setup(
    name='modelrun_series',
    version='0.1',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'Click',

    ],
    entry_points='''
        [console_scripts]
        modelrun_series=ripcas_dflow.scripts.modelrun_series:cli
    '''
)
