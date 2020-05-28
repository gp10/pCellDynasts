# pCellDynasts
Analyzing proliferating Cell Dynamics in squamous epithelial tissues.

Computational methods used to unravel the paradigm of cell proliferation in squamous epithelial tissues. This set of methods allows to respond to how many proliferating cell types exist, their renewal dynamics, and how lineage fates are orchestrated to guarantee tissue homeostasis.

Check the original publication where this methodology was applied for context details:
  > Piedrafita G, Kostiou V, Wabik A, Colom B, Fernandez-Antoran D, Herms A, Murai K, Hall BA, Jones PH (2020) A single-progenitor model as the unifying paradigm of epidermal and esophageal epithelial maintenance in mice. _Nat. Commun_ 11, 1429. https://doi.org/10.1038/s41467-020-15258-0

Tools described here can potentially be extended to disentangle the proliferating cell hierarchies and the mode of cell renewal in other tissues.

### Graphical abstract
![GraphicalAbstract](https://github.com/gp10/pCellDynasts/blob/master/Graphical_abstract_pCellDynasts.png)

### Overview
This code combines computational methods used for cell-fate model simulations, inference and fitting of experimental data, organized in three categories found as separate folders:
- **H2BGFP-ModalityTests**: inference on cell proliferation heterogeneity through modality tests. Statistical analysis of multimodality of H2BGFP dilution patterns across individual cells: hints on division rate homogeneity and thus on how many proliferating cell types exist.
- **H2BGFP-tcc-Inference**: inference on cell cycle time and division rate of proliferating cell populations. Cell division timing properties of the single-progenitor population are resolved from H2BGFP intensity distributions.
- **Clonal-Lineage-Inference**: inference on clonal dynamics. Division and differentiation fate probabilities as well as other kinetic properties are inferred from cell lineage tracing.

Documentation (Readme files) on scripts and code used for each specific purpose can be found in the corresponding folder.

### Requirements:
- Matlab (majority of code)
- Python
- R (for modality tests)

### Notes:
Figure and page references are now updated to match those in the publication in _Nat. Commun_ specified above. Note these do not match the order in the _bioRxiv_ preprint.
