
from dolfin import Function

def _Function_set(function, i, val):
    function.vector()[i] = val

def _Function_array(function):
    return function.vector().array()

Function.__setitem__ = _Function_set

Function.array = _Function_array