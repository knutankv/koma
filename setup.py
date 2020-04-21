import setuptools

with open("README.rst", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="koma-knutankv",
    version="0.1",
    author="Knut Andreas Kvåle",
    author_email="knut.a.kvale@gmail.com",
    description="Knut's Operational Modal Analysis Toolbox for Python. Matlab version available under './+koma/'.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/knutankv/koma",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)