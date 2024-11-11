from setuptools import setup, find_packages

setup(
    name='BSF',
    version='0.1',
    description='Fast light simulation using a beam-spread-function approach',
    author='David Berling',
    author_email='berling@ksvi.mff.cuni.cz',
    packages=find_packages(),
    install_requires=[
        'numpy', 
        'pandas',
        'scipy',
        'pyyaml',
        'matplotlib'
    ],
    test_suite='pytest',
    tests_require=['pytest'],
)

