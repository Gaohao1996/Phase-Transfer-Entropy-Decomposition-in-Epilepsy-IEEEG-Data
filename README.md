# Phase Transfer Entropy (TE) Decomposition in Epilepsy IEEEG Data

Functions for calculating TE Decomposition in IEEG data to find the redundant and synergistic paths beween certain source and target variables. The dictionary 
`PLV Functional Network` is used to visualize the functional network of the data based on their PLV values and do some analyses. The dataset used here is from [Statistical Analysis of Network Data](https://math.bu.edu/people/kolaczyk/datasets.html). Some of data has been preprocessed ( these can be identified in the `.mat` file which are), like removing the noise power frequency and the noise component from ICA ( only for seizure 1).  

## Transfer Entropy Decomposition
The main script `EZ_TEdecom` contains：
1) Oscillation  frequency band  filtered from data, i.e., Delta, Theta, Alpha, Beta (FIR filters are used here).
2) Time-frequecy power analyses to discover the potential EZ zone.
3) TE decompostion on Phase data. For the parameter of model order and delay in TE calculation, AIC is used to select the best order and delay.
   
Note: When testing TE with different source, target and conditioned variables, the covariance matrix could be ill conditioned at some time. In this case, try to use shrink covariance instead in function `cmi_ggg` to solve this problem. If it sitll doesn't work, you need to reduce the dimension of these variables.



## PLV Functional Network
The main script `PLV_and_functional_network` contains:
1) Calculate the PLV value between channels in certain frequency band (Delat, Theta, Alpha, Beta).
2) Signicance test of PLV value via P value and false correction to create functional networks.
3) Visualiztion and statistical analyses of functional networks in different time windows.

Note: This script aims to find the frequency band that is most related to the seizure activity so that we can descide to do TE decomposition in which band.
