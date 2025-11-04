from setuptools import setup, find_packages

setup(
    name="gaia",  # The name of the package
    version="0.1.0",  # The version of the package
    packages=find_packages(),  # Automatically discover and include all packages in the project
    install_requires=[  # External dependencies required for the project
        "numpy>=1.24.3",  # Required version of numpy
        "pandas>=1.3.5",  # Required version of pandas
        "scikit-learn>=1.2.2",  # Required version of scikit-learn
        "scanpy>=1.9.2",  # Required version of scanpy
        "geomstats>=2.5.0",  # Required version of geomstats
    ],
    python_requires=">=3.8",  # Python version requirement
    author="Jinpu Cai",  # Your name
    author_email="jinpucai99@gmail.com",  # Your email address
    description="GAIA: Geometric Analysis for single-cell RNA-seq",  # A short description of the package
    long_description="""GAIA provides tools for geometric analysis and visualization
                        of single-cell RNA-seq data using spherical embeddings and PCA.""",  # A detailed description of the package
    long_description_content_type="text/plain",  # The format of the long description (plain text)
    url="https://github.com/Carroll105/GAIA",  # URL to the project homepage or GitHub repository
    classifiers=[  # Classifiers to describe the project (useful for the package index)
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    keywords="single-cell RNA-seq, geometric analysis, PCA, scanpy, geomstats",  # Keywords related to the project for searchability
)
