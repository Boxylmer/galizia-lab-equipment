# # if we run a compiled C object, this will not be the method by which we enter
# from julia import Pkg
# Pkg.activate(".\\JuliaFunctions")
# from julia import JuliaFunctions
# JuliaFunctions.f(1, 2)
# print("done with julia")

import ctypes
from ctypes.util import find_library
from ctypes import *

path = os.path.dirname(os.path.realpath(__file__)) + '\\JuliaFunctions.dll'

# _lib = cdll.LoadLibrary(ctypes.util.find_library(path))
# hllDll = ctypes.WinDLL(path, winmode=0)
with os.add_dll_directory(os.path.dirname(os.path.realpath(__file__))):
    _lib = ctypes.CDLL(path, winmode=0)
print(path)
