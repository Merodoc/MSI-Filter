library("Cardinal")
path2 = paste("S043_Processed", ".imzML")
Experiment = readMSIDAta(path2)

register(SerialParam())
# Create a simulated Spectrum 
set.seed(2020)
mse <- simulateImage(preset = 1, npeaks=10, nruns=2,baseline=1)
mse

# Pixel Information stored in a Position Data Frame"""
pixelData(mse)

# extracts the pixel coordinate data frame
coord(mse)

# Vector of experimental runs
run(mse)[1:10]

#Feature Data
#m/z values 
featureData(mse)

#specifically the m/z vectors
mz(mse)[1:10]

#Image Data

#extract entries

iData(mse, "intensity")

#spectra() accessed first data matrix

spectra(mse)[1:5, 1:5]

#row of these matrices correspond to mass features
#columns to pixels

