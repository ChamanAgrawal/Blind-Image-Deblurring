# Software Requirements
MATLAB >= 2019a
CVX Package

# Code

## What is what
There are three main files (MATLAB LIVE SCRIPTS):
- ImageDeblurComparision.mlx : Live script for Image deblurring experiment.
- phase_transition.mlx : Live script for Phase Transition experiment.
- Robustness_comparision.mlx : Live script for Robustness experiment.

Other than these live scripts, there are four helper matlab scripts for ImageDeblurComparision.mlx
- DeblurAlgorithm1.m : Convex BD
- DeblurAlgorithm2.m : PRIDA
- DeblurAlgorithm3.m : Graph BD
- DeblurAlgorithm4.m : Non-convex BD

Then, there are four source code folders containing the actual implementation of the algorithms:
- Convex_src
- Non_convex_src
- Prida_src
- Graph_src

## How to run
Simply use the three MATLAB live scripts available to get all the presented results.

## Code Bibliography
- Non_convex_src     : Self-implemented
- DeblurAlgorithm1.m : Official implementation (https://github.com/CACTuS-AI/Blind-Deconvolution-using-Convex-Programming)
- DeblurAlgorithm2.m : Self-implemented
- DeblurAlgorithm3.m : GitHub implementation (https://github.com/BYchao100/Graph-Based-Blind-Image-Deblurring)
- DeblurAlgorithm4.m : Self-implemented
- ImageDeblurComparision.mlx : Mostly Self-implemented with some initial part from the mentioned Convex Official implementation.
- phase_transition.mlx       : Mostly from the mentioned Convex Official implementation with around 30% Self-implemented.
- Robustness_comparision.mlx : Mostly from the mentioned Convex Official implementation with around 30% Self-implemented.
- Non_convex_src     : Self-implemented
- Convex_src         : Official implementation (https://github.com/CACTuS-AI/Blind-Deconvolution-using-Convex-Programming)
- Prida_src          : Official implementation (https://github.com/sravi-uwmadison/prida)
- Graph_src          : GitHub implementation (https://github.com/BYchao100/Graph-Based-Blind-Image-Deblurring)




# Plots
All the plots used in the paper are stored in the 'plot' folder


