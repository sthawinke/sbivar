
# 0.1.0

- Package creation and first upload to Github

# 0.1.2

- Enforce valid names on data matrices with giveValidNames()

# 0.1.3

- Make the package more verbose, printing progress on fitting and
  testing

# 0.1.4

- Enforce x-y column names for the coordinates in sbivarMulti too

# 0.2.0

- Implement bivariate Moran’s I for single images
- Disallow sample matching and correlation for disjoint coordinate sets
- Scale by maximum values for bivariate Moran’s I

# 0.2.1

 - Modify weighting parameters for Moran's I
 - Add pseudocount for GAMs with gamma distribution
 - Stricter row and column name checking
 
# 0.99.1

- Preparing submission to Bioconductor

# 0.99.2

 - When Gamma fit fails, attempt negative binomial fit

# 0.99.3

 - Include Gaussian random field in GAM fitting to account for stochastic neighbour resemblance

# 0.99.4

 - Fix bugs in testGAM
 - Update Github readme

# 0.99.5

 - Fix bugs in plotGAM
 - Add feature names as columns in results

# 0.99.6

 - Switch to common names for metabolites in Vicari data
 - Switch off GP smoother in GAMs as default
 - In plotGAMsTopResults, use provided families as default
 
# 0.99.7

 - Rearrange matrix multiplications for GAMs to save memory
 - C++ speed-ups for GPs and Moran's I
 - Fixing some potential issues in fitLinModels and extractResultsMulti
 - Add Makevars for C++ Windows installation

# 0.99.8

 - Use bpworkers(bpparam()) to retreive number of registered cores rahter than bpparam()$workers
 
# 0.99.9

 - Call gtsat:variogram and gstat::vgm explicitly
 
# 0.99.10

 - Switch to SnowParam for multithreading on windows
 
# 0.99.10
 
 - No more multithreading in vignette building on windows
 
# 0.99.11

 - Add pseudocount for gamma family in plotGAM
 
# 0.99.12

 - Add option to scale point sizes to library sizes in case of normalization or offset modelling
 
# 0.99.13

 - Fix naming bug in extractResultsMulti
 
# 0.99.14

 - Keep finite correlation estimates 
 
# 0.99.15

 - Scale contributions to correlation by maximum absolute value in plotGAMs
 
# 0.99.16

 - Provide featuresX and featuresY to sbivarMulti
 
# 0.99.17

- Added r-universe github build
- Switch to tempfile for illustrating writeSbivarToXlsx()

# 0.99.18

 - Switch to lapply for buildGams
 
# 0.99.19

 - Estimate variograms only for features to be tested
 
# 0.99.20

 - Bug fixes when only one feature pair is tested
 
# 0.99.21

 - Round to significant digits in writeSbivarToXlsx()
 
# 0.99.22

 - Extend troubleshooting section in vignette
 - scaleBySampleSums default FALSE in plotPairSingle
 - New function for plotting the spatial distribution of the spot-wise sums
 
# 0.99.23
 
 - Bug in plotSpotSums
 
# 0.99.24

 - Citation file and reference to preprint
 
# 0.99.25

 - Remove unused offsets argument
 
# 0.99.26

 - Remove GP smooth, it is too concurve with the other smooth after all
 - Add GP error structure with mgcv::gamm instead, as an option

# 0.99.27

 - Switch to anndataR package rather than anndata and zellkonverter for reading in h5ad files
 - Fix minor bug in plotGAMs
 
# 0.99.28

  - More flexibility and better defaults for residual correlation structures
