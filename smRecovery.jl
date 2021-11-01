using MPIReco, JLD, PyPlot

####################
# load system matrix
####################
# Dataset store handling
# datadir = "./sm"
# store = MDFDatasetStore(datadir)

# system matrix
datadirSF2 = joinpath(calibdir(store), "2.mdf")
fSF2 = MPIFile(datadirSF2)
nx,ny,nz = calibSize(fSF2)
S2 = reshape( getSystemMatrix(fSF2,bgCorrection=true), nx, ny, :)
snr2 = MPIReco.getSNRAllFrequencies(fSF2)

# select a slice and apply SNR thresholding
snr_idx2 = findall(x->x>3.0, vec(snr2))[2:end]
S2 = S2[:,:,snr_idx2]

########################
# undersampled data
########################
rf = 6
numSamp = div(nx*ny,6)
#poisson disk
idx_pd = load("./samplingPatterns/pd_r$(rf).jld","idx")[1:numSamp]
# optimized patterns
idx_opt= load("./samplingPatterns/samplingIdxOpt.jld","idx")[1:numSamp]

# perform undersampling
y_pd = zeros(ComplexF64, length(idx_pd),size(S2,3))
y_opt = zeros(ComplexF64, length(idx_opt),size(S2,3))
for k=1:size(S2,3)
  y_pd[:,k] .= vec(S2[:,:,k])[idx_pd]
  y_opt[:,k] .= vec(S2[:,:,k])[idx_opt]
end

############
# SMRecovery
############
params = Dict{Symbol,Any}()
params[:shape] = (nx,ny)
params[:reg1] = "L1"
params[:sparseTrafo] = DCTOp(ComplexF64,(nx,ny),2)
params[:ρ_l1] = 0.3
params[:reg2] = "Nothing"
params[:iterationsInner] = 50
params[:iterations] = 10
params[:relTol] = 1.e-3
params[:absTol] = 1.e-4

params[:λ_l1] = 0.4^6
S2_pd = reshape( smRecovery(y_pd,idx_pd,params), nx, ny, :)

params[:λ_l1] = 0.4^5
S2_opt = reshape( smRecovery(y_opt,idx_opt,params), nx, ny, :)

################################
# plot some frequency components
################################
freqs=[24,35,38]
figure("recovered SMs", figsize=(6,6))
subplot(3,3,1); imshow(abs.(S2[:,:,freqs[1]])); title("k=$(snr_idx2[freqs[1]])"); ylabel("ref"); xticks([]); yticks([])
subplot(3,3,2); imshow(abs.(S2[:,:,freqs[2]])); title("k=$(snr_idx2[freqs[2]])"); xticks([]); yticks([])
subplot(3,3,3); imshow(abs.(S2[:,:,freqs[3]])); title("k=$(snr_idx2[freqs[3]])"); xticks([]); yticks([])
subplot(3,3,4); imshow(abs.(S2_pd[:,:,freqs[1]])); ylabel("PD"); xticks([]); yticks([])
subplot(3,3,5); imshow(abs.(S2_pd[:,:,freqs[2]])); xticks([]); yticks([])
subplot(3,3,6); imshow(abs.(S2_pd[:,:,freqs[3]])); xticks([]); yticks([])
subplot(3,3,7); imshow(abs.(S2_opt[:,:,freqs[1]])); ylabel("Opt"); xticks([]); yticks([])
subplot(3,3,8); imshow(abs.(S2_opt[:,:,freqs[2]])); xticks([]); yticks([])
subplot(3,3,9); imshow(abs.(S2_opt[:,:,freqs[3]])); xticks([]); yticks([])
tight_layout()
