from setuptools import setup
import versioneer

setup(
    name="ProbeDesign",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Design probes",
    author="Greg Gydush",
    author_email="ggydush@broadinstitute.org",
    packages=["probe_designer"],
    install_requires=[
        "numpy>=1.16.2",
        "pandas>=0.23.4",
        "click>=7.0",
        "biopython>=1.73"
    ],
    entry_points={
        "console_scripts": [
            "probe_designer = probe_designer.probe_designer:probe_designer"
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]
)
