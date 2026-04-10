from setuptools import setup, find_packages

with open("README.md", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="diffpeak",
    version="1.0.0",
    description="Differential peak analysis for CUT&Tag / ChIP-seq",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Hengbin Gao",
    author_email="hengbin.gao@uwa.edu.au",
    url="https://github.com/hengbingao/diffpeak",
    python_requires=">=3.8",
    packages=find_packages(),
    install_requires=[],
    entry_points={
        "console_scripts": [
            "diffpeak=diffpeak.analysis:run",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
    ],
    keywords="ChIP-seq CUT&Tag differential peaks bioinformatics",
)
