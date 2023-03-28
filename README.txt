README:

About:

MATLAB codes to read and process the crystal structure file (.xyz format) and calculate the inter/intra-octahedral distortion parameters (Octahedral bond length distortion (Δd), octahedral bond angle variance (𝛔2), and octahedral tilting and rotation angle (𝛳c and 𝛳ab),) in halide perovskites.

Requirements:

You need to have MATLAB software installed in your system.

Installation:

First copy the structure file and the Matlab code into the same folder. Then run the .m file with MATLAB. 

Uses:

We have provided the step-by-step procedure in the comment lines of each code. The usage of each available code is summarised below.

(1)	octahedral_elongation_and_compression.m: This code reads the crystal structure and calculates the B-X bond elongation and compression in each octahedron. To symmetrize the distortion data, one of the high symmetry axes is fixed along the z-direction. For this purpose, we have used the "rotation_mat_3D.m" code which rotates each octahedron by a given angle (defined by th, fi, and si in the code). 

(2)	distortion_parameters.m: This code reads the crystal structure and calculates the distortion parameters. The parameters (Δd, σ^2, 𝛳ab, and 𝛳c) for each structure are calculated for all inequivalent octahedra using the formula given in Eq. 2, 3, and 4 of the main text.


N.B. 

For further queries please contact nandab@iitm.ac.in and/or mayank77338@gmail.com.

Reference:
--------------------------


