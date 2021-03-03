# **Galizia Lab Equipment Code**

## Overview
This project can be broken up into categories for each available apparatus within the lab. Sections generally have tools used in previous studies and are developed as they are needed (i.e., usually not ahead of time). 

If you need a tool developed, a bug fixed, or feature added, please feel free to contact whoever maintains this at the time and they should be able to work with you. 


#Installation Notes
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
