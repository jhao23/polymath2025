# Ramsey Number R(J_5, K_5) Project Code

This repository contains the code for during the Polymath Jr program in 2025.

Generally, gluing problems (stored as compressed and pickled files) are created in the "pointed" python programs. For the case of degree 14, an extra "problemmaker" program is also provided as there are many gluing problem files that are created in this case.
These gluing problem files can then be read by the "gluing" programs which produce actual Ramsey graphs in (J_5, K_5) that are formatted in the graph6 format.
Then, the vertexExtension.cpp program attempts to add vertices to these graphs to see if any can be extended past the current lower bound. 

## Note: 
Running the "pointed" and "gluing" python programs will produce many hundreds of gigabytes of gluing problem/graph files on some of the bigger cases (as computing Ramsey numbers tends to be very computationally intensive).
I recommend that you set aside a lot of memory and break the computation into portions if you are to run the programs.

## Dependencies: 
The programs will depend on the external packages networkx, pynauty, and parsl. The vertex extension code is also compiled with C++20.

The reference graphs were obtained from Robert Fidytek's dataset of Ramsey graphs and code from Ash Kiel's graph enumeration algorithm (https://github.com/AshKila/J4K5Enumeration).

An alternative method for enumerating these graphs using a SAT solver is found here under Niccolo Turillo's github (https://github.com/nturillo/polymath2025gluing and https://github.com/nturillo/satjointsolver).
