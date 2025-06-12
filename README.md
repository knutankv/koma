![KOMA logo](https://raw.githubusercontent.com/knutankv/koma/master/koma-logo.svg)

What is KOMA?
=======================
KOMA is a package for operational modal analysis in Python. For additional details about the implementation of the covariance-driven stochastic subspace identification algorithm please refer to [5]. For automatic OMA and clustering analysis, please refer to [6]. More information and functionality will be added after publication of the cited paper.


Installation 
========================
Either install via PyPI as follows:

```
pip install koma-python
```

or install directly from github:

```
pip install git+https://www.github.com/knutankv/koma.git@master
```


Quick start
=======================
Import the relevant package modules, exemplified for the `oma` module, as follows:
    
```python
from koma import oma
```

For details, please refer to the examples. For code reference visit [knutankv.github.io/koma](https://knutankv.github.io/koma/).


![Shear frame](https://raw.githubusercontent.com/knutankv/koma/master/mode2.gif)


Examples
=======================
Examples are provided as Jupyter Notebooks in the [examples folder](https://github.com/knutankv/koma/tree/master/examples).

References
=======================
<a id="1">[1]</a> 
L HERMANS and H VAN DER AUWERAER. MODAL TESTING AND ANALYSIS OF STRUCTURES UNDER OPERATIONAL CONDITIONS: INDUSTRIAL APPLICATIONS. Mechanical Systems and Signal Processing, 13(2):193–216, mar 1999. URL: http://www.sciencedirect.com/science/article/pii/S0888327098912110, doi:http://dx.doi.org/10.1006/mssp.1998.1211.

<a id="2">[2]</a>
Peter Van Overschee and Bart De Moor. Subspace identification for linear systems: theory, implementation, applications. Kluwer Academic Publishers, Boston/London/Dordrecht, 1996.

<a id="3">[3]</a>
Carlo Rainieri and Giovanni Fabbrocino. Operational Modal Analysis of Civil Engineering Structures. Springer, New York, 2014.

<a id="4">[4]</a> Brad A. Pridham and John C. Wilson. A study of damping errors in correlation-driven stochastic realizations using short data sets. Probabilistic Engineering Mechanics, 18(1):61–77, jan 2003. URL: http://www.sciencedirect.com/science/article/pii/S0266892002000425, doi:10.1016/S0266-8920(02)00042-5.

<a id="5">[5]</a> Knut Andreas Kvåle, Ole Øiseth, and Anders Rønnquist. Operational modal analysis of an end-supported pontoon bridge. Engineering Structures, 148:410–423, oct 2017. URL: http://www.sciencedirect.com/science/article/pii/S0141029616307805, doi:10.1016/j.engstruct.2017.06.069.

<a id="6">[6]</a> K.A. Kvåle and Ole Øiseth. Automated operational modal analysis of an end-supported pontoon bridge using covariance-driven stochastic subspace identification and a density-based hierarchical clustering algorithm. IABMAS Conference, 2020.

Citation
=======================
Zenodo research entry:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10277112.svg)](https://doi.org/10.5281/zenodo.10277112)

Support
=======================
Please [open an issue](https://github.com/knutankv/koma/issues/new) for support.

