import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="apscale_blast", # Replace with your own username
    version="2.0.5",
    author="Till-Hendrik Macher",
    author_email="macher@uni-trier.de",
    description="Advanced Pipeline for Simple yet Comprehensive AnaLysEs of DNA metabarcoding data - BLAST application",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://pypi.org/project/apscale_blast/",
    packages=setuptools.find_packages(),
    license='MIT',
    install_requires=[
                    'Bio >= 1.7.1',
                    'biopython >= 1.84',
                    'joblib >= 1.4.2',
                    'ete3 >= 3.1.3',
                    'numpy',
                    'pandas >= 2.2.2',
                    'pyarrow >= 16.1.0',
                    'xmltodict >= 0.14.2',
                    'streamlit_extras'
                      ],
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.10',
    entry_points={
        "console_scripts": [
            "apscale_blast = apscale_blast.__main__:main",
        ]
    },
)

