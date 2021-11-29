using Pkg
#Pkg.add("EllipsisNotation")
#Pkg.add("BenchmarkTools")
#Pkg.add("StaticArrays")
#Pkg.add("CPUTime")
Pkg.add("Plots")
Pkg.add("GR")
Pkg.add("PyPlot")
Pkg.add("PyCall")
Pkg.add("HDF5")
Pkg.add("Test")
Pkg.add("LaTeXStrings")
Pkg.add("FFTW")

ENV["PYTHON"] = "/usr/bin/python3"
Pkg.build("PyCall")

Pkg.precompile()
