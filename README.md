# LES-TT
Large Eddy Simulation for Transcritical Turbulence

Transcritical turbulent flows are governed by the compressible Navier-Stokes equations along with a real-gas equation of state. Their computation is strongly susceptible to numerical instabilities and requires kinetic-energy- and pressure-equilibrium-preserving schemes to yield stable and non-dissipative scale-resolving simulations. Building upon a recently developed kinetic-energy- and pressure-equilibrium-preserving discretization framework based on transporting a pressure equation, the objectives of this work are to (i) derive a filtered set of equations suitable for large-eddy simulation, and (ii) characterize the properties of the resulting subfilter-scale terms by performing a priori analyses of transcritical wall-bounded turbulence direct numerical simulation data. The filtering operation leads to three unconventional subfilter-scale terms that emerge from the pressure equation and require dedicated modeling. 
The subfilter-scale stress tensor is dissected in terms of magnitude, shape and orientation based on an eigendecomposition analysis, and compared with existing subfilter-scale models. 
Therefore, the scripts embedded in this repository aim to carry out the these a propri analysis.

Based on filtered DNS the a priori analysis allow the user to perform a term-by-term analysis of each contributor from momentum and pressure transport equations and also the non-linear equation of state.
The filtered data can be evaluated in terms of:
· Wall-normal ensemble-average metrics such as turbulent kinetic energy, velocity, temperatures, first- and second-order statistics, instanteneous contourplots.
· Resolved (filtered) terms of the Navier-Stokes equations and unresolved (Sub-filter scale (SFS))
· Eigen decomposition of stress tensor
· Characterization of state-of-art SFS models for stress tensor such as: Classical Smagorinsky, Dynamic Smagorinsky, Anisotropic minimum-dissipation (AMD), WALE, Sigma and Similarity models
· Eigen decomposition of the SFS stress tensor models
· Characterization of stress tensor trace models: Yoshizawa and Vreman
· Correlation coefficient between filtered DNS and SFS models
· Evaluation of Taylor expansion and ILA models for the non-linear equation of state
· Ploting functionalities for either option

The proper execution of these scripts require the download of HPCFS Solver as it uses part of common utilities, the Peng-Robinson equation of state and high-pressure coefficient models, convective and viscous terms functionalities.

Multiple DNS snapshots are also provided to performed the filtered DNS analysis. Instead, the data container which contains the filtered data based on the average of these various datasets is also provided. However, due to the limitation of large files in GitHub, the data is stored locally on clusters at Universitat Politècnica de Catalunya · BarcelonaTech (UPC), and will be provided upon making proper arrangements with the requesters.

The MATLAB directory requires a Data and Figures folder.

1) DNS_filter.m

This script (i) filters the DNS snapshots based on the filter type and filter width selected and (ii) performs a ensemble-average based on multiple DNS snapshots.

In addition a csv file is created with the wall-normal (XZ) ensemble average quantities and a mat file with the filtered fields for each filter width so that it can be used for future post-processing.

Finally the base plots are also populated if bPlot boolean is true (Plot_filter.m) when the ensemble-average turbulent kinetic energy, 1st and 2nd order statistics and snapshots are plot.
  
2) Unclosed_terms.m

Given a selected filter width it loads the *.mat file and assesses each term-by-term analysis for momentum, pressure and equation of state.

It also populates the wall-normal ensemble-average plots (Plot_UnclosedTerms.m).

Finally, if the user runs the script "EigenDecomposition.m" the EigenDecomposition of the filtered DNS data loaded will be computed including the PDF histogram plots for trace, barycentric maps for the shape and polar maps for the orientation.


3) SFS_Models.m

This function requires a previous run of Unclosed_terms.m so that the results of a given filter width are loaded in the Workspace.

Hence, given these filter width results, the various SFS Models are exectuted. The user needs to be consistent on the filter selected, in this case the baseline is the CLDF_Box.

This post-process calculates the features for each model, the FG-level quantities and additional filtered required to perform trace models, stress tensor models and non-linear equation of state models.

Trace models:
- Yoshizawa
- Vremann

Stress tensor models:
- Classical Smagorinsky
- Dynamic Smagorinsky
- WALE
- Sigma
- AMD
- Similarity

Equation of state models:
- Taylor expansion
- ILA


This function also performs the Eigendecomposition of the SFS models (EigenDecomposition_SFS_Model.m) and generates the plots accordingly.

4) SFS_closure_expression.m

This code evaluates some of the clousure expressions for the unresolved terms which reported significant weight with respect to their resolved counterpart. Hence, smilarity-based models are evaluated and proposes.

