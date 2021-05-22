# Blind-Image-Deblurring
Blind Deconvolution is the problem of recovering
two unknown signals `w` and `x` from their circular convolution
with or without noise. This problem arises in many areas
such as astronomy, wireless communication, etc.  in the form of 
Image processing, signal processing, and so on. One of the major challenges in the
Blind Deconvolution problem is its non-convex nature. It has
been observed that this problem suffers from a large number of
local minimas when posed as an optimization problem.
In this project we have implemented four different algorithms to
solve this problem and have tried to compare their performance
under various settings. Two of the algorithms are studied over
general signal vectors while the other two are specifically studied for the
Image Deblurring problem (non-negative signals).

# Software Requirements
- MATLAB >= 2019a
- CVX Package

# Code overview
There are three main files (MATLAB LIVE SCRIPTS):
- ImageDeblurComparision.mlx : Live script for Image deblurring experiment.
- phase_transition.mlx : Live script for Phase Transition experiment.
- Robustness_comparision.mlx : Live script for Robustness experiment.

Other than these live scripts, there are four helper matlab scripts and four source code folders:
- DeblurAlgorithm1.m (Convex BD)
- DeblurAlgorithm2.m (PRIDA)
- DeblurAlgorithm3.m (Graph BD)
- DeblurAlgorithm4.m (Non-convex BD)
- Convex_src
- Non_convex_src
- Prida_src
- Graph_src

# How to run
Simply use the three MATLAB live scripts available to get all the presented results.

# Code Bibliography
- Non_convex_src     : Self-implemented
- Convex_src         : [Official implementation](https://github.com/CACTuS-AI/Blind-Deconvolution-using-Convex-Programming)
- Prida_src          : [Official implementation](https://github.com/sravi-uwmadison/prida)
- Graph_src          : [GitHub implementation](https://github.com/BYchao100/Graph-Based-Blind-Image-Deblurring)
