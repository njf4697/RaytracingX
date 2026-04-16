#RaytracingX

A relativistic raytracing thorn in the CarpetX framework.

Written by Noelle Feist <njf4697@rit.edu> with code adapted from Jhon Sebastián Moreno Triana <github: Jh0mpis> in collaboration with Manuela Campanelli, Yosef Zlochower, Liwei Ji, and Scott Noble.

This thorn creates a grid of mass-less particles (sometimes referred to as photons in the code) set to move in different directions originating at a single point, then integrates their corresponding geodesics backwards until they either reach the simulation boundaries, an event horizon, or the photosphere (WIP). The geodesic integration is modified from and has [GInX](https://github.com/Jh0mpis/GInX) as a prerequisite. That point is then taken to be the starting point for integration of the relativistic radiative transfer equation, and the specific intensity is integrated back to the camera (WIP).


This thorn creates initial position and velocities for geodesic integration based on camera parameters defined by the user. The implementation follows https://arxiv.org/pdf/1410.7775.
