This is the directory containing strategy files that we designed to project 2 of the Programming & Problem Solving course http://www.cs.columbia.edu/~kar/4444f16/

Notes are kept here: https://docs.google.com/document/d/1eVZND3tshZ2UHWBc0euW4Xj16scnhxJJJy2aCAIugQ8/edit

The files themselves are not meaningful unless they are compiled together with the simulator, for those who have the simulator and want to change/add files under this directory, please keep in mind we should rename this "project2-g4" to "g4" under your local directory ./slather/

Tips to help you run/compile the code:

To run or compile the simulator, cd into the folder above slather

To (re)compile the simulator on Unix & Mac OS X:   

    $ javac slather/sim/*.java

To (re)compile the simulator on Windows:          

    $ javac slather\sim\*.java

To run the simulator:  

    $ java slather.sim.Simulator <arguments>
    
The simulator arguments are:

    -g <group name(s) separated by spaces, e.g. g1 g2 g3>
    -t <integer number of turns pherome lasts, default 10>
    -d <double vision radius of cell, default 2>
    -s <side length of board, default 100>
    --fps <integer specifying fps for gui display>
    --gui
    --verbose

*Note that due to the hardcode in the simulator, -g <> must be the last parameter! Otherwise you cannot run (details of dealing with arguments could be seen at simulator.java)

For example, we should usually run like this (our g4 vs two dummy g0): 

    $ java slather.sim.Simulator --gui -g g0 g0 g4
