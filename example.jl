##################################
# Load Packages and Activate Env #
##################################
using Pkg 
Pkg.activate(".")
Pkg.instantiate()

using MPIReco, JLD, PyPlot

###############
# download data
###############
datadir = "./sm" #artifact"data"
store = MDFDatasetStore(datadir)

#####################
# run example scripts
#####################
include("smSampling.jl")
include("smRecovery.jl")