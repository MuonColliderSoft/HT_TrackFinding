# convertBIB

The ROOT macro convert_BIBRecoHits.C runs on the tracker's RecoHit collections stored in an LCTuple, adjusts the hit positions to lay on the simplified tracker geometry used in the HT_toy, and dumps the converted hits into a flat text file.


# HT_toy

The HT toy consists of two executables:

HTArrayTraining: training of the Hough-transform array. 

EventGeneration: generation of tracks with the BIB overlaid. 


# plotTracks

The ROOT macro PlotTracks.cpp reads the track candidates from the file PlotData.txt and displays the candidate hits. 
