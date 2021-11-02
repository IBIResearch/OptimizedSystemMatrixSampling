using MPIReco, JLD, SensorSelection, PyPlot

# estimate measurement noise from the training data
function varMeas(x::Matrix{T},rs::Int64) where T
  roi = vec(x[1:rs,1:rs].-mean(x[1:rs,1:rs]))
  append!(roi,x[end-rs+1:end,1:rs].-mean(x[end-rs+1:end,1:rs]))
  append!(roi,x[1:rs,end-rs+1:end].-mean(x[1:rs,end-rs+1:end]))
  append!(roi,x[end-rs+1:end,end-rs+1:end].-mean(x[end-rs+1:end,end-rs+1:end]))

  return sum(abs.(vec(roi)).^2)/(length(roi)-1)
end

####################
# load system matrix
####################
# Dataset store handling
# uncomment the next two lines if running this script separately
# datadir = artifact"MDFStore"
# store = MDFDatasetStore(datadir)

# system matrix
datadirSF1 = joinpath(calibdir(store), "1.mdf")
fSF1 = MPIFile(datadirSF1)
nx,ny,nz = calibSize(fSF1)
S1 = reshape( getSystemMatrix(fSF1,bgCorrection=true), nx, ny, :)
snr1 = MPIReco.getSNRAllFrequencies(fSF1)

# select a slice and apply SNR thresholding
snr_idx1 = findall(x->x>3.0, vec(snr1))[2:end]
S1 = S1[:,:,snr_idx1]

# normalize the frequency components
for k=1:size(S1,3)
  S1[:,:,k] ./= norm(S1[:,:,k])
end

# chose frequency components for training
Random.seed!(1234)
p = shuffle(collect(4:2:size(S1,3)))
St = S1[:,:,p[1:30]]  # training data

########################################################
# setup up measurement operators and covariance matrices
########################################################

# mean amplitudes of measurements
st_amp = zeros(nx,ny,size(St,3))
for k=1:size(St,3)
  st_amp[:,:,k] .= sqrt.( 1/size(St,3)*dropdims(sum(abs.(St[:,:,:]).^2,dims=3),dims=3) )
end
st_amp = reshape(st_amp,nx*ny,:)

# measurement matrix
dctOp = DCTOp(ComplexF64,(nx,ny))
H = Matrix(adjoint(dctOp))
Ht = collect(transpose(real.(H))) # build transpose of measurement matrix

# DCT representation of training data
St_dct = zeros(ComplexF64,nx*ny,size(St,3));
for k=1:size(St,3)
    St_dct[:,k] .= dctOp*vec(St[:,:,k])
end

# obtain sparse measurement matrix
thresh = 1.e-1  # portion of the total number of sparse coefficients to use for optimization 
numCoeff = floor(Int64,thresh*nx*ny)
Ht_sub = zeros(numCoeff,size(Ht,2),size(St,3))
for k=1:size(St,3)
  Sk_dct = vec(abs.(St_dct[:,k]))
  idx_sig = maxIndices(Sk_dct,numCoeff)
  Ht_sub[:,:,k] .= Ht[idx_sig,:]
end

# estimate measurement noise for each frequency component
measNoise = [varMeas(St[:,:,k],3) for k=1:size(St,3)]
r = zeros(nx*ny,size(St,3))
for k=1:size(St,3)
  r[:,k] .= fill(measNoise[k],nx*ny)
end

######################
# Forward optimization
######################
λ=3.0                             # proportionality factor for the prior covariance matrices
numSamples = floor(Int64, nx*ny)  # number of samples to add to the sampling pattern
w = zeros(Bool,nx*ny)             # initial sampling pattern (empty)
m = SensorSelection.Experiment(Ht_sub, r, Ht_sub, λ.*st_amp.^2,w)
@time idx_opt = optSamplingFW!(m,numSamples,4,180)

#############################
# plot some sampling patterns
#############################
# precomputed Poisson disk patterns
msk2_pd = zeros(nx,ny); msk2_pd[load(datadir*"/samplingPatterns/pd_r2.jld","idx")] .= 1
msk6_pd = zeros(nx,ny); msk6_pd[load(datadir*"/samplingPatterns/pd_r6.jld","idx")] .= 1
msk10_pd = zeros(nx,ny); msk10_pd[load(datadir*"/samplingPatterns/pd_r10.jld","idx")] .= 1
# optimized patterns
msk2_opt = zeros(nx,ny); msk2_opt[idx_opt[1:div(nx*ny,2)]] .= 1
msk6_opt = zeros(nx,ny); msk6_opt[idx_opt[1:div(nx*ny,6)]] .= 1
msk10_opt = zeros(nx,ny); msk10_opt[idx_opt[1:div(nx*ny,10)]] .= 1

figure("sampling Patterns", figsize=(8,6))
subplot(2,3,1); imshow(msk2_pd); title("r=2"); ylabel("PD"); xticks([]); yticks([])
subplot(2,3,2); imshow(msk6_pd); title("r=6"); xticks([]); yticks([])
subplot(2,3,3); imshow(msk10_pd); title("r=10"); xticks([]); yticks([])
subplot(2,3,4); imshow(msk2_opt); ylabel("Opt"); xticks([]); yticks([])
subplot(2,3,5); imshow(msk6_opt); xticks([]); yticks([])
subplot(2,3,6); imshow(msk10_opt); xticks([]); yticks([])
tight_layout()