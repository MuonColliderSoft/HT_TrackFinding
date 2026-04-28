

# Reconstruction.cpp

## Main Features

- Reads configurable event samples from input track files.
- Optionally overlays Beam-Induced Background (BIB) hits.
- Uses a Hough Transform Array for pattern recognition.
- Selects candidate cells in parameter space.
- Fits candidate tracks.
- Applies chi-square quality cuts.
- Suppresses duplicate candidates.
- Measures reconstruction efficiency and candidate multiplicity.
- Produces ROOT histograms for offline analysis.

## Physics / Detector Quantities Reconstructed

The program reconstructs track parameters such as:

- `phi`
- `eta`
- `1/pT`
- `z0`
- `t0`
- `beta`

These are compared to generated truth values during debugging and validation.

## Workflow

### 1. Initialize Hough Transform Array

At startup the program:

- creates the HTA object
- initializes histograms
- checks array dimensions
- loads training / mapping data from file

### 2. Read Events

For each event:

- reads one or more signal tracks
- adds optional BIB background hits
- builds a combined hit list

### 3. Fill Pattern Recognition Array

Seed hits are inserted into the HTA structure using the selected fill mode.

### 4. Find Best Cells / Candidates

The code identifies:

- best populated cells
- candidate cells in parameter space

### 5. Fit Tracks

Each candidate is fitted and assigned:

- chi-square
- reduced chi-square
- fitted track parameters

Good fits passing quality cuts are accepted.

### 6. Remove Duplicates

Equivalent solutions mapping to the same final cell are marked so only one survives.

### 7. Store Histograms

The program writes many ROOT histograms, including:

#### Event / Reconstruction

- number of candidates
- number of found tracks
- number of fits

#### Fit Quality

- total chi-square
- reduced chi-square
- good-fit chi-square
- layers used in fit

#### Signal Hit Distributions

- X
- Y
- Z
- T
- XY occupancy
- Z vs radius

#### Background Hit Distributions

- X
- Y
- Z
- T
- XY occupancy
- Z vs radius



# CutOptimization.cpp

The program learns sensor-by-sensor cut windows from simulated events, then measures how many background hits survive those cuts. 

## Main Goals

1. Read simulated events containing detector hits.
2. Group hits by sensor.
3. Build distributions for three observables:
   - X1
   - X2
   - T (time)
4. For each sensor and coordinate:
   - determine min / max
   - compute quantile limits
   - derive optimized low/high cuts
5. Save diagnostic histograms:
   - PDF files for quick inspection
   - ROOT histograms for later browsing/editing
6. Read an external background-hit file.
7. Apply the learned cuts sensor-by-sensor.
8. Produce detailed summary tables of background rejection.

## Input Data

The code uses standard Muon Collider software objects such as:

- Track
- Hit
- Event
- TrackReader
- BibFileReader

It reads:

### Signal / reconstruction sample
Used to build the coordinate distributions and derive cuts.

### Background sample
Used afterward to test the cuts.

## Per-Sensor Optimization

For each sensor with sufficient statistics, the code analyzes:

- X1 distribution
- X2 distribution
- Time distribution

From these, it stores:

- Minimum
- Maximum
- Lower quantile
- Upper quantile
- Final CutLow
- Final CutHigh

These values are organized in a map keyed by sensor name and coordinate, for example:

- Sensor_X1
- Sensor_X2
- Sensor_T

## Histograms Produced

For every sensor and coordinate, the program can create:

- PDF plots in the `hists/` folder
- ROOT histograms inside the output ROOT file

Typical names:

- Sy..._X1
- Sy..._X2
- Sy..._T

This allows later inspection with ROOT TBrowser.

## Background Rejection Study

After optimization, the code reads a background-hit sample and counts, for each sensor:

- Nin   : total hits entering this sensor
- NX1   : hits passing X1 cut
- NX2   : hits passing X2 cut
- NT    : hits passing time cut
- NAll  : hits passing all three cuts simultaneously

A final table is printed with totals across all sensors.

## Statistical Prediction Section

The code also estimates expected triple-cut acceptance assuming independence:

P(pass all) = P(X1) × P(X2) × P(T)

and compares:

- Expected accepted hits
- Actual accepted hits

This helps diagnose correlations among variables.

## Typical Workflow

1. Run reconstruction sample.
2. Build optimized cuts automatically.
3. Inspect histograms.
4. Apply cuts to background sample.
5. Compare expected vs observed rejection.

## Output Summary

The program prints:

### Cut Table

For every sensor and coordinate:

- Min
- Max
- QMin
- QMax
- CutLow
- CutHigh

### Background Table

For every sensor:

- Nin
- NX1
- NX2
- NT
- NAll

### Prediction Table

For every sensor:

- pX1
- pX2
- pT
- Expected N
- Actual N

## Dependencies

- C++17
- ROOT
- Muon Collider reconstruction software classes
- Local helper modules included by source

