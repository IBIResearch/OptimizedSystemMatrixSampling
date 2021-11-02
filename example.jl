##################################
# Load Packages and Activate Env #
##################################
using Pkg 
Pkg.activate(".")
Pkg.instantiate()

using MPIReco, JLD, PyPlot, LazyArtifacts

###############
# download data
###############
datadir = artifact"MDFStore"
store = MDFDatasetStore(datadir)

#####################
# run example scripts
#####################
include("smSampling.jl")
include("smRecovery.jl")