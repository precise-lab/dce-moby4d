import os
import dolfin as dl


source_directory =  os.path.dirname(os.path.abspath(__file__) )
fname = "source.cpp"

with open(os.path.join(source_directory, fname), "r") as cpp_file:
        cpp_code  = cpp_file.read()
        include_dirs = [".", source_directory]
        cpp_module = dl.compile_cpp_code(cpp_code, include_dirs=include_dirs)

SlitBeam = cpp_module.SlitBeam
ConeBeam = cpp_module.ConeBeam
BoundaryQ0Source = cpp_module.BoundaryQ0Source