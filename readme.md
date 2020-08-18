# vtk_tools 

This package provides useful utilities for working with vtk output files. It's current focus is on handling tasks related to cell ordering and underlying shape functions for different cell types. 

## installation 

none yet.... 

## requirements 

While `setup.py` (doesn't exist yet) contains requirements, it is important to stress that this package requires vtk version of `9.0.1` or later. At present, the latest vtk version on conda is 8.2, so if you already have installed vtk from conda (or from source), you will need to upgrade. The easiest way is to use pip: 

`pip install vtk`

If you already have a vtk version installed, you can try `pip install vtk -upgrade` but this will likely fail if the previous version of vtk was not installed with pip, in which case the previous version must be uninstalled frist. 

Also, note that while this package requires vtk `9.0.1` or later, it can process output that was created with lower vtk versions.
