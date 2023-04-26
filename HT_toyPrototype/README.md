# HT_toy

The HT toy package consists of three executables:

- *ConvertBIBHits*: it runs on the tracker's RecoHit collections stored in an LCTuple, adjusts the hit positions to lay on the simplified tracker geometry used in the HT_toy, and saves the converted hits into a flat text file.
- *HTArrayTraining*: training of the Hough-Transform array.
- *EventGeneration*: generation of tracks with the BIB overlaid.

Instructions:

- to compile the executables:\
  ```make```

- to convert the BIB tracker hits:\
  ```./ConvertBIBHits <input LCTuple file> <output text file>```

- to run the Hough-Transform array training:\
  ```./run_HTtraining.sh```

- to generate events:\
  ```./EventGeneration```
  

# plotTracks

The ROOT macro *PlotTracks.cpp* reads the track candidates from the file *PlotData.txt* and displays the candidate hits. 

Usage:\
```root -l PlotTracks.cpp```
