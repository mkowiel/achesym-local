# Achesym-local

Achesym - standardized placement of macromolecular models in the unit cell, command line tool

This is a simple script that allows to use the achesym program on local machine.

# Citing

    Kowiel, M., Jaskolski, M. & Dauter, Z. (2014). 
    ACHESYM: an algorithm and server for standardized placement of macromolecular models in the unit cell. 
    Acta Cryst. D70, doi:10.1107/S1399004714024572. 

# Server
 
To use Achesym web server please visit: http://achesym.ibch.poznan.pl/

# Installation 

1. Install Python. The script is compatible with Python 2.7.
2. Install cctbx (https://github.com/cctbx/cctbx_project)
3. `cd achesym`
4. Execute script `cctbx.python achesym.py -analyze input.pdb output.pdb`

# Usage

    achesym.py input_file.pdb output_file.pdb [structure_factors.cif output_structure_factors.cif]
                     [-analyze] [-pack] [-group "A,B;C,D"]
    
    -analyze: Find the most compact assembly (with "common volume" analysis)
    -pack: Find the most compact assembly (with bounding boxes) - deprecated
    -group "A,B;C,D": force chain grouping. In the example chains A and B will be
                      treat like one object and chains C and D as a second object
                      group separator ";"
                      chain separator ","
    
    Program requires cctbx, and can be run by cctbx aware python intepreter
    For example:
    #> cctbx.python achesym.py -analyze in.pdb out.pdb 

