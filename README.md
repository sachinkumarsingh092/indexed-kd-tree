# Quad-Hash Matching

The main objective of project is to create quad-hashes from the five brightest stars in the reference and query cataloge.

To make proper hashes, we use a healpix like grid `struct` to include all the stars with their respective brightness. These stars are then used to create quad-hashes (which consists of 4 stars) and then made a local coordinate system and the hashes are calculated from them. 

These hashes are then stored in a kd-tree for quick retrievals and nearest neighbour calculations.

This repository contains just the files that are to be finally mearged in GNU-Astro.