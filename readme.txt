README                      10/09/2008
http://code.google.com/p/diffusion-mri
------------------------------------------

------------------------------------------
| Requirements for Using the Python Code |
------------------------------------------

    Platforms:
        The Python code should work on all platforms which are supported by Python, NumPy and SciPy.

    Requirements:
    (1) Python (version >=2.4 )
        Most versions of Linux and Mac OS already have python preinstalled. Windows users can download
        Python from http://www.python.org for free. Please make sure that Python is put in your path,
        so that typing 'python' in the system command line starts the Python interpreter.
    (2) NumPy (version >=1.0) and SciPy (version >=0.6.0)
        In addition to the standard python installations, both NumPy and Scipy have to be installed so
        that neither typing 'import numpy' nor 'import scipy' in the python interpreter raise an error.
        For how to obtain and install NumPy/SciPy, please see http://www.scipy.org/


-------------------------------
| How to Test the Python Code |
-------------------------------

    Data preparation:
        To use this program, you need to have the input diffusion MRI data,
        acquisition information ready and prepare a configuration file (.ini) to specify path and various
        parameters. The details of the input data and configuration file will be covered in the rest of
        this README.

    Example usage with the provided sample data and configuration file:
      (1) Assume the program you checked out from svn chunk is in c:\diffusion-mri.
      (2) To test the Python code, please go to the data subdirectory in the system command line and type
          'python ../Python/reconstruction.py example.ini mow'. The first argument gives the configuration
          file and the second argument specifies the reconstruction method. 'mow' stands for mixture of
          Wisharts,  other two available methods are 'qbi' and 'dot'.
      (3) You can also test the Python code in the Python interpreter or an interactive Python shell. if
          IPython (http://ipython.scipy.org/) is installed, which is strongly recommended, please try
          the following sequence of commands in IPython:
             cd c:\diffusion-mri\Python
             from reconstruction import MOWReconstructor
             cd ..\data
             reconstructor = MOWReconstructor()
             reconstructor.set_config_file('example.ini')
             reconstructor.reconstruct()


---------------------------------------------------------
| About Configuration File and Input/Output Data Format |
---------------------------------------------------------

    Configuration file:
        We use the standard ini file to store the information about where to read input data, acquisition
        parameters, user-specified model parameters for different reconstruction methods. Please refer to
        the example.ini in the data subdirectory which is self-explained.

    Data format:
        Both the input diffusion image data and the reconstructed data are stored in the mhd file format.
        For input diffusion image data, the last dimension is the number of diffusion gradients. For the
        reconstructed data, the last dimension is the number of spherical harmonics coefficients.

    Other files:
        Other parameters like diffusion gradients, b-matrices, spherical harmonics basis matrices
        and tessellation schemes are stored in plain text files like 81vectors.txt, 321vectors.txt


-----------------------------------------------------------------
| Last modified: 10/09/2008                                     |
| If you have any questions, please contact bing.jian@gmail.com |
-----------------------------------------------------------------
