from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='bsf_light',
    version='1.0',
    description='Beams Simple and Fast LIGHT simulation (bsf_light) - Fast light simulation for brain tissue using a beam-spread-function approach',
    long_description=long_description,
    long_description_content_type="text/markdown",
    author='David Berling',
    author_email='berling@ksvi.mff.cuni.cz',
    packages=find_packages(),
    install_requires=[
        'numpy', 
        'scipy',
        'pyyaml',
        'matplotlib'
    ],
    test_suite='pytest',
    tests_require=['pytest'],
)
