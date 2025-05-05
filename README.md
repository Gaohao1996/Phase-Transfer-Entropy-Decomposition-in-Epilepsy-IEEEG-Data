# Phase Transfer Entropy Decomposition in Epilepsy IEEEG Data

Functions for calculating Transfer Entropy Decomposition in IEEG data with the popurse of finding the redundant and synergistic path beween source and target variables. Also the dictionary PLV Functional Network is used to visualize the functional network of the data based on their PLV values and do some analyses. The dataset used here is from . The data has been preprocessed, like removing the noise power frequency and the noise component from ICA ( only for seizure 1).  

## Transfer Entropy Decomposition
The main script `EZ_TEdecom` contains：
1) Oscillation  frequency band  filtered from data, i.e., Delta, Theta, Alpha, Beta (FIR filters are used here).
2) Time-frequecy power analyses to discover the potential EZ zone
3) TE decompostion on Phase data



## PLV Functional Network
The main script  contains:
1) Calculate the PLV value between channels in certain frequency band (Delat, Theta, Alpha, Beta).
2) Signicance test of PLV value via P value and false correction to create functional networks
3) Visualiztion and statistical analyses of functional networks in different time windows 
