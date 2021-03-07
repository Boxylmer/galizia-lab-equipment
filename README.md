# **Galizia Lab Equipment Code Repository**

## Overview
This project can be broken up into categories for each available apparatus within the lab. Sections generally have tools used in previous studies and are developed as they are needed (i.e., usually not ahead of time). 


## Installation Notes

The python package `llvmlite` doesn't currently support python `python 3.9`, so you need to install a virtual environment with `python 3.8` or below.

For `python 3.8`, the easiest way to do this is the `py` launcher
-  Using the `py` launcher, in contrast to using `python`, will actually interpret your command and understand what versions of python you have installed, rather than just loading whatever shows up in your #path

1.  Install `python 3.8` (do not add it to your path)


2. Create the virtual environment
    
    `py -3.8 -m venv env`


3. Activate the environment 
    
    `env\Scripts\activate`


4. install the project dependencies
    
    `pip install requirements.txt`

## Features and to-do list

#### Core code
- [x] Create symbolic analytical error propagation equation generator
- [x] Create Monte-Carlo error propagation generator
- [ ] Create comparison and visualization tools for generators

#### Vapor Sorption Apparatus 
- [x] Create symbolic vapor sorption apparatus representation
- [x] Write function to export symbolic vapor sorption apparatus to a python function
- [ ] ~~Test and / or determine if exported functions can be optimized in Numba~~ **Not easily**
- [x] Calculate a single monte carlo step from given input data
- [ ] Calculate an entire isotherm using the Monte-Carlo method
- [ ] Write test methods for analytical representation and numerical representation comparisons
  
   
   - (long term goals)
- [ ] Define a format for generating and reading template files to quickly analyze sorption isotherm data

#### Gas Sorption Apparatus
- [ ] Create python function representation
- [ ] Figure out the best way to import the PREoS code sustainably
   - [ ] Check how easily this can simply be implemented in *Julia*?

#### Permeation Cell Apparatus
- [ ] Create python function representation
- [ ] Write test methods for MCEP
