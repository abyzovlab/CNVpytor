from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name = 'CNVpytor',
    version = '1.0a1',
    packages = ['cnvpytor'],
    package_dir={'cnvpytor': 'cnvpytor'},
    package_data={'cnvpytor': ['data/*.pytor']},
    description="Python extension of CNVnator",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/abyzovlab/CNVpytor",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"
    ],
    entry_points = {
        'console_scripts': [
            'cnvpytor = cnvpytor.__main__:main'
        ]
    })
