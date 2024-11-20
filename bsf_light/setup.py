from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='bsf_light',
    version='1.0.1',
    description='Beams Simple and Fast LIGHT simulation (bsf_light) - Fast light simulation for brain tissue using a beam-spread-function approach',
    long_description=long_description,
    long_description_content_type="text/markdown",
    author='David Berling',
    author_email='berling@ksvi.mff.cuni.cz',
    packages=find_packages(),
    python_requires=">=3.13",
    install_requires=[
        'numpy>=2.1.3', 
        'scipy>=1.14.1',
        'pyyaml>=6.0.2',
        'matplotlib>=3.9.2'
    ],
    test_suite='pytest',
    tests_require=['pytest'],
)
