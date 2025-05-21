# Phase Transfer Entropy (TE) Decomposition in Epilepsy IEEEG Data

Functions for calculating TE Decomposition in IEEG data to find the redundant and synergistic paths beween certain source and target variables. The dictionary 
`PLV Functional Network` is used to visualize the functional network of the data based on their PLV values and do some analyses. The dataset used here is from [Statistical Analysis of Network Data](https://math.bu.edu/people/kolaczyk/datasets.html). Some of the data have been preprocessed through removing the power frequency noise and the noise component from ICA ( only do ICA for seizure 1), which can be identified by the name of the variable in `.mat` files ( clean or raw data).  

## Transfer Entropy Decomposition
The main script `EZ_TEdecom` contains：
1) Oscillation filtered from data via FIR filters, i.e., Delta, Theta, Alpha, Beta.
2) Time-frequecy power analyses (Short-time fourier tranfermation, STFT) to discover the potential EZ zone.
3) TE decompostion on Phase data. For setting the parameter of model order and delay in TE calculation, Akaike Information Criterion (AIC) is used here.
   
Note: When testing TE with different source, target and conditioned variables, the covariance matrix could be ill conditioned at some time. In this case, you should try to use shrink covariance instead in the function `cmi_ggg` to solve this problem. If it sitll doesn't work, you need to reduce the dimension of these variables.



## PLV Functional Network
The main script `PLV_and_functional_network` contains:
1) Calculate the PLV value between channels in certain frequency band (Delat, Theta, Alpha, Beta).
2) Signicance test of PLV value via P value and false correction (`fdr_bh`) to create functional networks.
3) Visualiztion and statistical analyses of functional networks in different time windows.

Note: This script aims to find the frequency band that is most related to the seizure activity so that we can descide to do TE decomposition in which band.
