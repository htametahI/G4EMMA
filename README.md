# S2223 Simulations

## Table of Contents

- [Technologies](#Technologies)
- [Run](#run)
- [Input Files](#input-files)
- [Output Files](#output-files)
- [Notes](#Notes)





## Technologies
Simulations are built with: 
- MacOs Tahoe 26.2 (M1 Chip)
- `ROOT` 6.36.06
- `GEANT4` 11.3.2
- `CMake` 4.2.1
- Apple `clang` version 17.0.0 (clang-1700.6.3.2)
- `f2c` static library available at `f2c/lib/libf2c.a`

## How to use 
- Clone the repo: https://github.com/htametahI/G4EMMA.git
- `cd G4EMMA`
- simply call `./cmake_script.sh`. You will have to modify the paths becasue they are hardcoded. 



## Folder Structure

- `src/`, `include/`: main simulation source code
- `UserDir/UserInput/`: run input files (`beam.dat`, `reaction.dat`, `centralTrajectory.dat`, etc.)
- `UserDir/Results/`: simulation outputs (ROOT + text data)
- `macros/`: Geant4 macro files
- `root_macros/`: ROOT analysis scripts/macros


## Run

From the repository root:

```bash
./EMMAapp
```


## Input Files

At startup, the application reads:

- `UserDir/UserInput/beam.dat`
- `UserDir/UserInput/reaction.dat`
- `UserDir/UserInput/centralTrajectory.dat`

## Output Files

Right now I have configured the Simulation to output .dat data files. Each row of the data represents a particular event, with different data fields specified by columns. This enables easy visulization in a Python Notebook. 

e.g. `S3_triton_ring_observables.dat`

The script `check_results.ipynb` implements data plotting and fitting. 

## Notes

In the simulation I am currently killing all tracks before EMMA to save computing time. To disable this, comment out`theTrack->SetTrackStatus(fStopAndKill)` in `EMMASteppingAaction.cc`

Currently the main output data file is: `S3_triton_ring_observables.dat`
In `Exc.csv`, I've selected a few states for EMMA transmission efficiency calculations. 


