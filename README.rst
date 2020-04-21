.. image:: koma-logo.svg

KOMA
=======================

KOMA is a package for operational modal analysis, with core functionality available both in Python and MATLAB languages. For additional details about the implementation of the covariance-driven stochastic subspace identification algorithm please refer to [1]. Data-SSI and FDD are implemented in the MATLAB version of KOMA. For automatic OMA and clustering analysis, please refer to [7]. More information and functionality will be added after publication of the cited paper. 

Please refer to the `documentation <https://gist.github.com/1855764>`_ for more details.


Installation and usage
-----------------------

Python
......................

Either download the repository to your computer and install, e.g. by **pip**

.. code-block::

   pip install .


or install directly from the python package index.

.. code-block::

   pip install koma


Thereafter, import the package modules as follows:
    
.. code-block:: python

    from koma import oma, modal, clustering


MATLAB
..............
Download or clone repository. Folder containing package root has to be in the MATLAB path:

.. code-block:: matlab

   addpath('C:/Users/knutankv/git-repos/koma/');

Ideally this is done permanently, such that it is always accessible. Then, the package can be
imported with the following command:

.. code-block:: matlab

   import koma.*

Now all the subroutines of the package are accessible through

.. code-block:: matlab

    koma.function_name

E.g., to use the function `covssi.m` located at .../+koma/+oma/ the following syntax is applied:

.. code-block:: matlab

   [lambda,phi,order] = koma.oma.covssi(data, fs, i, 'order', order);

Functions inside `private` folders are accessible only from functions
inside the folder at the root of the `private` folder.

Citation
=======================
Please cite both the original research paper for which the code was developed for and `the code <https://zenodo.org/record/NUMBER>`_ if used in research. 

Support
=======================
Please `open an issue <https://github.com/knutankv/koma/issues/new>`_ for support.

Acknowledgements
=======================
Jos van der Geest contribution to Mathworks Central File Exchange is used to produce error bars in stability plots of MATLAB function stabplot.m.

References
=======================
[1] Knut Andreas Kvåle, Ole Øiseth, and Anders Rønnquist. Operational modal analysis of an end-supported pontoon bridge. Engineering Structures, 148:410–423, oct 2017. URL: <http://www.sciencedirect.com/science/article/pii/S0141029616307805>, doi:10.1016/j.engstruct.2017.06.069.
