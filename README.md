# TimeResolvedNetworks

In case of using this code in your work, please cite “Spatiotemporally flexible subnetworks reveal the quasi-cyclic nature of integration and segregation in the human brain”, Strindberg et al 2021


TimeResolvedNets.m creates time-resolved spatiotemporally flexible subnetworks as described in Strindberg et al 2021 with associated time-courses. It also creates state vectors recurrence matrices and estimates flexibility/modularity of individual brain areas.


TimeResolvedNets.m can be run using either the empirical mode decomposition (EMD)-algorithm to derive intrinsic mode functions (IMFs) for subsequent calculation of subnetworks. For this choice use “band = 0”. The user will subsequently be prompted to chose which IMF to use. If instead a band pass filter is preferred, chose “band = 1”. The user will then be promoted to set the lower and upper limit of the filter.


TimeResolvedNets.m have two subroutines: SNCs.m and PhaseIntegrationSNC.m


SNCs.m calculates the SNCs based on the seed pairs calculated in TimeResolvedNets.m. The user chose the desired size of final SNCs (n = 4, n = 6, n = 8, n = 10, n = 12). If one instead want to start with a set of predefined single areas as seeds, these can be entered in vector form (varargin). The SNCs will then have odd numbers (n = 3, n = 5, n = 7, n = 9, n = 11).


PhaseIntegrationSNC.m calculates the phase coherence within a subset of the SNCs as well as randomly composed components of the same size prior to and during integration (assignment into the same community of the areas of the SNC/random component). The user is prompted to chose the number of SNCs and random components that should be included in the calculation. 


All code except PhaseIntegrationSNC.m can be run for either a single subject or multiple subjects. PhaseIntegrationSNC.m is only adapted to run for a minimum number of two subjects.


Except from the standard MATLAB library the following code is also used and should be downloaded prior to use:  

The Louvain community detection algorithm from the Brain Connectivity Toolbox: https://sites.google.com/site/bctnet/measures/list

Permutation Test by Laurence R Krol: https://github.com/lrkrol/permutationTest

Mantel test by Enrico Glerean: https://figshare.com/articles/software/Mantel_test/1008724

The BrainNet viewer was used for visualization of networks: www.nitrc.org/projects/bnv (Xia et al., 2013). 

MaxVox.mat is used for the visualisation of subnetworks (SNs) and meta-networks (MNs).

AreaParcelA236_N9.mat contains the network number as defined by the parcellation (Schaefer200 7 network parcellation combined with the subcortical portion of the Brainnetome atlas) in the first column and area number in the second column.  Details on the parcellation can be found in Parcellation.xlsx. The parcellation is contained in Parcellation236mixed.zip.

NetNames.mat matches network number with network name. 

If a different parcellation is used all three mat-files needs to be modified accordingly for the full code to run properly. To update MaxVox.mat first update HCP236.node according to the desired parcellation. Then open MaxVox.mat in the BrainNetViewer and simply change to the new node file and save.

For questions and bug fixes please contact:  marikastrindberg@gmail.com or marika.strindberg@ki.se
