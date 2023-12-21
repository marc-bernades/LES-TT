# LES-TT
Large Eddy Simulation for Transcritical Turbulence

1) DNS_filter.m

This script (i) filters the DNS snapshots based on the filter type and filter width selected and (ii) performs a ensemble-average based on multiple DNS snapshots.

In addition a csv file is created with the wall-normal (XZ) ensemble average quantities and a mat file with the filtered fields for each filter width so that it can be used for future post-processing.

Finally the base plots are also populated if bPlot boolean is true (Plot_filter.m) when the ensemble-average turbulent kinetic energy, 1st and 2nd order statistics and snapshots are plot.
  
2) Unclosed_terms.m

Given a selected filter width it loads the *.mat file and assesses each term-by-term analysis for momentum, pressure and equation of state.

It also populates the wall-normal ensemble-average plots (Plot_UnclosedTerms.m)


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

  

