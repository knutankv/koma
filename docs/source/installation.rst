Installation and usage
=======================

Python
-----------
Either download the repository to your computer and install, e.g. by **pip**

.. code-block::

   pip install .


or install directly from the python package index.

.. code-block::

   pip install git+https://www.github.com/knutankv/koma.git@master
   

Thereafter, import the package modules, exemplified for the `oma´ module, as follows:

.. code-block:: python

    import koma.oma


MATLAB
-----------
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

E.g., to use the function :file:`covssi.m` located at .../+koma/+oma/ the following syntax is applied:

.. code-block:: matlab

   [lambda,phi,order] = koma.oma.covssi(data, fs, i, 'order', order);

Functions inside `private` folders are accessible only from functions
inside the folder at the root of the `private` folder.
