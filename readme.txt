Prerequisites
-------------

(1) Python (version 2.4 or higher)
(2) Numpy (version 1.0 or higher)
(3) Scipy (version 0.6.0)


Example Usage
--------------
python reconstruction.py example.ini mow
python reconstruction.py example.ini dot
python reconstruction.py example.ini qbi

Please see the comments in the python source files and the example.ini
for more details.



Display the Reconstruction Results
------------------
The visualization tool vis_spharm.sav is written in Interactive Data Language (IDL) and 
uses a graphical user interface to config the parameters. It has been tested with IDL 
version 6.1 and higher, and it runs under the IDL Virtual Machine (IDL VM) 


- IDL VM should be installed first;

	 The IDL Virtual Machine is a free runtime version of IDL available
from ITT Visual Information Solutions at http://www.ittvis.com/idlvm/

- Windows users:
       1. Double click vis_spharm.sav;
       2. Use the "Choose File" button to select file to display (_real.flt file);
       3. Edit box "Visualization Dimension" can be used to change size;
       4. Use "Visualize" button to start visualization;

- Linux users:
       1. In the directory of file vis_spharm.sav, use command:
               $ idl -vm=vis_spharm.sav
       2. The same as 2,3,4 steps for Windows users;

- Hint:
       1. _real.flt and _imag.flt file should be in the same directory.
       2. vis_spharm.sav can only display data in 2D lattice. 3D volume data needs to be sliced first.





