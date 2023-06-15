Talha Hussain, MyMeshFairingModelingProjectAssignmentFour, Computer Science
Computer Graphics Project Assignment 4 Mesh Fairing Modeling

Websites:
To complete my Computer Graphics Project Assignment 4 Mesh Fairing Modeling, I received helpful inspiration from the following websites:
Implicit Fairing of Irregular Meshes using Diffusion and Curvature (Caltech Joural Article) - http://multires.caltech.edu/pubs/ImplicitFairing.pdf
For implementing the BLAS functions needed for the Linear Solver. https://netlib.org/blas/  

Usage:
A Makefile has been included with the following commands:
make smoothing - Runs the compilation using the g++ compiler

make clean - cleans the intermediate files and compiled target

After running the Makefile command "make smoothing", choose one of the following make commands below to run and then make clean. 
Each test set other than test-assignment-image will run with four strength and iteraction count combinations:
    strength 1 iterations 1
    strength 3 iterations 1
    strength 1 iterations 50
    strength 3 iterations 50

make test-assignment-image - Runs the smoothing on bunny.obj in an attempt to generate the same output of the bunny depicted in the assignment documentation for
Computer Graphics Project Assignment 4 Mesh Fairing Modeling.
make test-basic - Runs the four strength and iteration combinations without using any flags.
make test-cotangent - Runs the four strength and iteration combinations using cotangent weights flag.
make test-biharmonic - Runs the four strength and iteration combinations using biharmonic/bilaplacian flag.
make test-preserve - Runs the four strength and iteration combinations using the preserve volume flag.
make test-implicit - Runs the four strength and iteration combinations using the implicit linear solver with an eps of 0.001.
make test-loop - Runs the four strength and iteration combinations using the loop subdivide flag with a subdivision count of 1.
make test-all - Runs all test commands one after the other.

Summary:
Mesh fairing is a complex topic with a simple result. After working on the Computer Graphics Project Assignment 4 Mesh Fairing Modeling, I learned that Mesh
Fairing is a mesh which has a smoother appearance than originally inputted. This is accomplished by slowly reducing the acuteness of the angles between vertices
causing the mesh to not be as angular.

I started with the basic Laplacian Calculation as described in the assignment by accumulating the differences between vertices to create an overall direction of
travel and a weight indicating how far to travel along that vector. This produced a robust result with enough iterations and strength to produce a very smooth
looking mesh.

For the Contangent Weights Extra Credit Part of Computer Graphics Project Assignment 4 Mesh Fairing Modeling, I calculated the weight not by how far the points were
but by the ratio of how long an edge was compared to its neighbor edges. This allows for a more nuanced weighting system in exchange for a much longer runtime
because the weights need to be recalculated for every iteration.

Next, I implemented the volume preservation as described in the attached Caltech journal article "Implicit Fairing of Irregular Meshes using Diffusion and Curvature"
by using the technique of calculating the triangle normal and distance from the center point of the object to determine the internal volume then multiplying each
point by the ratio of size change an iteration had produced.

Implicit integration was added using the Linear Solver code and implicit.cpp example given. I modified the mult function to use the basic Laplacian calculator to
perform the smoothing and calculating the eps error until the iteration was complete.

Finally loop subdivision was accomplished by performing a loop over the triangles and analyzing the creation of mid edge points before rebuilding the triangles
using the previous and new point data. Although the process of Loop Subdivision was a bit slow, it however resulted in a correct subdivided triangle and a good
quality of resulting mesh.

In conclusion with that done and all tests producing results, I was glad with what I compiled for my Computer Graphics Project Assignment 4 Mesh Fairing Modeling.
