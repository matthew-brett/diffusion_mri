README                      09/27/2008
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


---------------------------------------
| Requirements for Using the C++ Code |
---------------------------------------

    Platforms:
        This C++ program has been tested on Windows XP, SUSE Linux 10.0.

    Build from the C++ source code:
      (1) CMake is used in the building process from the C++ source code. If CMake is not installed on 
          your system, please download it from: http://www.cmake.org and see the documentation there.
      (2) We use the vnl library in VXL for the numerical computation in this program. So you need to 
          first install the VXL from http://vxl.sourceforge.net/ also by CMake. VXL is a collection of 
          C++ libraries designed for computer vision research and implementation. However, for building 
          our program, only vnl and vcl are needed.


---------------------------------------
| How to Test the C++ and Python Code |
---------------------------------------

    Data preparation:
        To use this program (either Python or C++ code), you need to have the input diffusion MRI data, 
        acquisition information ready and prepare a configuration file (.ini) to specify path and various 
        parameters. The details of the input data and configuration file will be covered in the rest of 
        this README.

    Example usage with the provided sample data and configuration file:
      (1) Assume the program you checked out from svn chunk is in c:\diffusion-mri.  
      (2) To test the C++ program, first make sure the executable 'reconstruction' (or 'reconstruction.exe') 
          is successfully built from the source code,then go to the data subdirectory in the system command 
          line and run the command "reconstuction example.ini mow". The first argument gives the configuration 
          file and the second argument specifies the reconstruction method. 'mow' stands for mixture of 
          Wisharts,  other two available methods are 'qbi' and 'dot'.
      (3) To test the Python code, please go to the data subdirectory in the system command line and type 
          'python ../reconstruction.py example.ini mow'. Similarily, 'mow' can be replaced by  'qbi' or 'dot'.
      (4) You can also test the Python code in the Python interpreter or an interactive Python shell. if 
          IPython (http://ipython.scipy.org/) is installed, which is strongly recommended, please try 
          the following sequence of commands in IPython:
             cd c:\diffusion-mri\Python
             from reconstruction import MOWReconstructor
             cd ..\data
             reconstructor = MOWReconstructor('example.ini')
             reconstructor.reconstruct()


---------------------------------------------------------
| About Configuration File and Input/Output Data Format |
---------------------------------------------------------

    Configuration file:
        We use the standard ini file to store the information about where to read input data, acquisition
        parameters, user-specified model parameters for different reconstruction methods. Please refer to 
        the example.ini in the data subdirectory which is self-explained. 

    FLT data format:   
        Both the input diffusion image data and the reconstructed data are stored in the  binary FLT format.
        A flt file starts with a header structure
        +-------------------------------------------------------------------+
        | dim (int) | size (int[dim]) | data_type (int) | data_length (int) |
        +-------------------------------------------------------------------+
        followed the binary data section.

        Notes:
        (1) The number of dimensions is stored in the first 4 bytes as integer type.
        The range of dimension should be from 1 to 4. It is very likely that the byte
        order used by the file is big-endian rather than the little-endian on your PC.
        So swapping may be required to get the right number.
        (2) The current program only recognizes float data type which is encoded as 4.
        For complex data, the data_type field should be 8.

        Example:
        If 'data.flt' contains 81 diffusion MR images, each of which has size 128x128x64,
        then the header organzied in 4-bytes-unit should look like
             4,  128, 128, 64, 81, 4, 128*128*64*81
        and right after the header are 128*128*64*81 float numbers stored by stacking 81
        float-valued images of size 128*128*64.

    Other files:
        Other parameters like diffusion gradients, b-matrices, spherical harmonics basis matrices
        and tessellation schemes are stored in plain text files like 81vectors.txt, 321vectors.txt

-------------------------------------------
| Visualization of the Reconstructed Data |
-------------------------------------------

    vis_spharm.sav and its requirement:

        A visualization tool with graphical user interface 'vis_spharm.sav' is included in the 'data' 
        subdirectory. It was written by Evren Ozarslan in Interactive Data Language (IDL). To run this 
        program, the IDL Virtual Machine (version 6.0 or later) has to be installed first.  The IDL 
        Virtual Machine (IDLVM) is a free runtime version of IDL available from ITT Visual Information 
        Solutions at http://www.ittvis.com/idlvm/.

        The following instructions assume that IDLVM is already installed.


    How to use the vis_spharm.sav:
        (1) In Windows system, double click the file 'vis_spharm.sav'. In Linux, go to the directory 
            containing the file 'vis_spharm.sav', type the command 'idl -vm=vis_spharm.sav'. This will 
            pop-up a little GUI that allows you to specify the inputfile and visualization parameters.
        (2) Click the "Choose File" button to select the reconstruction result file (_real.flt file) to 
            be visualized. Note, for each _real.flt, there must be a corresponding _imag.flt file in 
            the same folder.
        (3) Click the "Visualize" button to start visualization        
        (4) Some words about the options that you can play with in this GUI:
            * The number "Visualization Dimension" is the size of square box for rendering a single voxel.
            * The number "Maximum order" is the highest spherical harmonics order that the program will 
              read from the input data
            * The background color can be changed by giving RGB values. For example, 255,255,255 gives 
              a white background color.
            * If the checkbox 'Overlay on GA' is checked, the generalized anisotropy (GA) index 
              (see Evren et al. MRM05) will be computed and encoded in color overlayed on each voxel.
            * If the checkbox 'Save PNG' is checked, the visualization result will be saved to a PNG 
              file given in the 'PNG File name' box.
        (6) Notes:
            * The checkbox 'Color' does not function properly, so please do not check the checkbox 'Color'.
            * Currently 'vis_spharm.sav' can only display data in 2D lattice. In order to visualize the
              reconstruction data from 3D volumes, the resulting flt files may have to be sliced into 
              2D data in order to be processed by the 'vis_spharm.sav' program.


-----------------------------------------------------------------
| Last modified: 09/27/2008                                     |
| If you have any questions, please contact bing.jian@gmail.com |
-----------------------------------------------------------------
