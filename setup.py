import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="APSCALE blast", # Replace with your own username
    version="1.0.0",
    author="Till-Hendrik Macher",
    author_email="macher@uni-trier.de",
    description="Advanced Pipeline for Simple yet Comprehensive AnaLysEs of DNA metabarcoding data - Graphical User Interface",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/apscale_blast",
    packages=setuptools.find_packages(),
    license='MIT',
    install_requires=['pySimpleGUI >= 4.15.2',
                      'apscale >= 1.5.5',
                      'demultiplexer >= 1.1.0',
                      'lastversion >= 1.3.3',
                      'boldigger >= 1.2.5',
                      'plotly >= 4.9.0',
                      'requests >= 2.25.1',
                      'requests_html >= 0.10.0',
                      'kaleido >= 0.0.3',
                      'xlsxwriter >= 3.0.2',
                      'tqdm >= 4.60.0',
                      'ete3 >= 3.1.2',
                      'biopython >= 1.78',
                      'xmltodict >= 0.12.0'
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
            "apscale_blast = apscale_gui.__main__:main",
        ]
    },
)
