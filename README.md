# Ramsey Number R(J_5, K_5) Project Code

This repository contains the code for during the Polymath Jr program in 2025.

Generally, gluing problems (stored as compressed and pickled files) are created in the "pointed" python programs. For the case of degree 14, an extra "problemmaker" program is also provided as there are many gluing problem files that are created in this case.
These gluing problem files can then be read by the "gluing" programs which produce actual Ramsey graphs in (J_5, K_5) that are formatted in the graph6 format.
Then, the vertexExtension.cpp program attempts to add vertices to these graphs to see if any can be extended past the current lower bound. 

Note: Running the "pointed" python programs will produce many hundreds of gigabytes of gluing problem files on some of the bigger cases (as computing Ramsey numbers tends to be very computationally intensive).
I recommend that you set aside a lot of memory or break the computation into portions if you are to run the programs.
