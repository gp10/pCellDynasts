# H2BGFP_tcc_Inference
Computational methods to estimate cell division rate & individual cell-cycle period distributions (tcc) from experimental H2BGFP dilution data.

### Main scripts
- **Analysis-H2BGFPdil_histograms_overview.py** : script used for plotting experimental keratinocyte H2BGFP distributions over time (the user can select the data set of interest and the script returns detailed histograms separated per individual mice or field-of-view).
- **Analysis_H2BGFPdil_theorModelPredictions.m** : script used for exploring theoretical H2BGFP decay distribution profiles under different models of cell behavior (SP model vs. SC-CP and 2xSC model paradigms).  
- **Analysis-H2BGFPdil-tccDist-inference.m** : main script to infer tcc, by comparing simulated vs. experimental H2BGFP dilution data (the user can select the data set of interest and the average division rate and best tcc distribution fits displayed).
- **Analysis-H2BGFPdil-tccDist-inference-Mascre2012.m** : script used for tcc parameter inference and SP model goodness-of-fit on cell division data from Mascre et al (2012).
- **Analysis-H2BGFPdil-tccDist-inference-Sada2016.m** : script used for tcc parameter inference and SP model goodness-of-fit on cell division data from Sada et al (2016).

### Dependencies
- Datasets : contains the experimental H2BGFP dilution data from the esophagus and different skin territories referred to in the manuscript.
- MonteCarloSimulator-SP-BasalCloneDynamics-CellProlif.m : function used to run time course simulations of basal-layer clone sizes and cell proliferation under the SP model.
- MonteCarloSimulator-SCCP-BasalCloneDynamics-CellProlif.m : function used to run time course simulations of basal-layer clone sizes and cell proliferation under the SC-CP model.
- MonteCarloSimulator-2xSC-TotalCloneDynamics-CellProlif.m : function used to run time course simulations of clone sizes and cell proliferation under the 2xSC model from Sada et al (2016).
- reshape_BackskinDataFormat.m : reformats the original back-skin data (coming from 3 different experiments) into a compact form similar to the one used for other data sets.
- ABCrejection-tccDist-inference.m : runs main ABC algorithm for tcc parameter estimation (based on experimental vs. simulation goodness-of-fit on H2BGFP distributions).
- ABCrejection-tccDist-inference-BasalCellProlif.m : variant of the previous ABC algorithm for tcc parameter estimation to accomodate it for contrasting distributions of the No. of basal cell divisions (e.g. Mascre 2012 data).
- ABCrejection-tccDist-inference-SuprabasalCellDiv.m : variant of the ABC algorithm for stratification/shedding waiting time parameter estimation when considering No. of cell divisions experienced by suprabasal (spinous) cells (Sada 2016 data).
- ABCrejection-set-TolThres-acceptRate.m : tests goodness-of-fit on H2BGFP distributions for random sets of parameter values to define an appropriate tolerance threshold for main ABC algorithm.
- ABCrejection-set-TolThres-acceptRate-BasalCellProlif.m : variant of the previous to test goodness-of-fit on distributions in the No. of cell divisions and define an appropriate tolerance threshold for variant ABC algorithm (e.g. on Mascre 2012 data).
- ABCrejection-set-TolThres-acceptRate-SuprabasalCellDiv.m : variant to test goodness-of-fit on distributions in the No. of cell divisions in suprabasal (spinous) compartment, so as to define appropriate tolerance threshold for variant ABC algorithm (Sada 2016 data).
- MonteCarloSimulator-SP-BasalCellProlif-H2BGFPdil.m : runs simulations of basal cell turnover and H2BGFP dilution under the SP model hypothesis.
- MonteCarloSimulator-SP-BasalCellProlif.m : runs simulations of basal cell proliferation under the SP model hypothesis (ignoring the H2BGFP dilution process).
- MonteCarloSimulator-SP-SuprabasalCellDiv.m : runs simulations of suprabasal cell replacement as basal cell proliferation proceeds (ignoring the H2BGFP dilution process).
- ODEs-2xSC-CellProlif-EXPtcc.m : ODEs describing cell proliferation and turnover rates (No. of division rounds in each compartment) under the 2xSC model hypothesis (assumption of exponential waiting time processes).
- ODEs-2xSC-CellProlif-GAMtcc.m : ODEs describing cell proliferation and turnover rates (No. of division rounds in each compartment) under the 2xSC model hypothesis (assumption of gamma (k=2) waiting time processes).
- SelectModelParamVal.m : function used to call to specific tcc parameter sets.

### Requirements
- Python 3.6.5	(Spyder with IPython 6.4.0 was used)
- Matlab R2016b
