# convertBIB

The ROOT macro convert_BIBRecoHits.C runs on the tracker's RecoHit collections stored in an LCTuple, adjusts the hit positions to lay on the simplified tracker geometry used in the HT_toy, and dumps the converted hits into a flat text file (BIBdata.txt).

Usage:
root -b -q convert_BIBRecoHits.C+\(\"BIB_LCTuple.root\"\) >& out_convertBIB.log


# HT_toy

The HT toy consists of two executables:
- HTArrayTraining: training of the Hough-Transform array.
- EventGeneration: generation of tracks with the BIB overlaid.

Instructions:
- to compile the executables:
  make

- to run the Hough-Transform array training:
  ./run_HTtraining.sh

- to generate events:
  ./EventGeneration
  

# plotTracks

The ROOT macro PlotTracks.cpp reads the track candidates from the file PlotData.txt and displays the candidate hits. 

Usage:
root -l ../plotTracks/PlotTracks.cpp
