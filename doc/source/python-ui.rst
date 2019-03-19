.. _python-ui:

Python UI
=========

The Python User Interface is a python module named :code:`pycompwa`, built with 
`pybind11 <https://pybind11.readthedocs.io/en/stable/index.html>`_.
It is the recommended way to use ComPWA, since the user benefits from the python ease of use.

.. note::
   Because the Python UI calls the c++ code in the background, you do not have to worry about speed.
   It runs just as fast ;)

The Python UI enables you to perform all of the tasks needed for your partial wave analysis: 

- load and create and intensity from a description
- load or generate data samples
- perform fits
- save & visualize results

On how to use the Python UI, please refer to the :ref:`examples section<examples>`.
Below you can find the code documentation of **pycompwa**.

.. automodule:: pycompwa
   :members:
   :undoc-members:
   :show-inheritance: