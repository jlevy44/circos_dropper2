from setuptools import setup

setup(
    name='circos_dropper',
    version='0.1',
    py_modules=['circos_dropper'],
    install_requires=[
        'numpy',
        'pandas',
        'pysam==0.10.0',
        'pybedtools',
        'pyfaidx',
        'pathos',
        'jcvi'
    ],
    entry_points="""
        [console_scripts]
        circos_dropper=circos_dropper:circosdrop"""
)