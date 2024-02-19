import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="koma-knutankv",
    version="1.2.3",
    author="Knut Andreas Kvåle",
    author_email="knut.a.kvale@gmail.com",
    description="Knut's Operational Modal Analysis Toolbox for Python.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/knutankv/koma",
    packages=setuptools.find_packages(),
    install_requires=['plotly', 'pandas', 'numpy', 'pyvista[jupyter]>=0.38.1',
                      'scikit-learn', 'trame', 'ipywidgets', 'pyvistaqt'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6'
)
