# HT_toyPrototype

The HT_toyPrototype package consists of four executables:

- *ConvertBIBHits*: it runs on the tracker's RecoHit collections stored in an LCTuple, adjusts the hit positions to lay on the simplified tracker geometry used in the HT_toy, and saves the converted hits into a flat text file.
- *HTArrayTraining*: training of the Hough-Transform array.
- *EventGeneration*: generation of tracks with the BIB overlaid.
- *PlotTracks*: it reads the track candidates from the file *PlotData.txt* and displays the candidate hits.

All the configurable parameters are set in the header file *Parameters.h*.


### Instructions:

- to compile the executables:\
  ```make```

- to convert the BIB tracker hits:\
  ```./ConvertBIBHits <input LCTuple file> <output text file>```

- to run the Hough-Transform array training:\
  ```source run_HTtraining.sh```

- to generate events:\
  ```./EventGeneration```

- to display the candidate hits:\
  ```./PlotTracks <input text file>```
