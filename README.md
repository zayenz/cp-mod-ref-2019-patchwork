# State Representation and Polyomino Placement for the Game Patchwork

This repository contains the code for the paper "State Representation
and Polyomino Placement for the Game Patchwork" by Mikael Zayenz
Lagerkvist. The paper will be presented at the 18th workshop on
Constraint Modelling and Reformulation co-located with the 25th
International Conference on Principles and Practice of Constraint
Programming 2019.

The code is located in the directory `code/` and is a CMake C++ 17
project.

## Building and running the experiments

[Gecode 6.2.0](https://www.gecode.org/ ) is required to be
installed on the system before building. The code has only been tested
on a Mac, it is likely that another way of discovering and linking
Gecode might be needed for other platforms, currently the CMake
`find_library` method is used for finding and setting up Gecode.

First, create a new build directory `code/build`. From within `build`, use 
```
$ cmake ..
$ make
``` 
to build the application.

There will be two executables produced,
`src/patchwork-test-placements` and
`src/patchwork-test-game-play-statistics`. The former runs the main
experiment, testing all placement strategies, while the latter runs
some simple game play experiments for game play statistics.


