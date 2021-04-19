This is a repository of code for the research titled "Applying DSA to the 2018-2020 Ebola Epidemic in Eastern DRC"

# Dynamic Survival Analysis Instructions
The data used for these demos are labled drc1_demo.csv, drc2_demo.csv, drc3_demo.csv.  This dataset consists of randomly chosen infection times, with recovery times created by adding a unifromly chosen number from uniform distribution 0 to 16.  This individual level infection and recovery time data can then be used in Dynamical Survival Analysis by running the following r scripts below.  The simulations are done using rstan.

## If recovery data is available ...
In the case that you have data for individual recovery time, as well as infection time, use the script and plot files beginning with "DSA".  Read in the data you would like to use, with the onset time first, followed by the column for recovery time. Note that the times for infection and recovery are continuous data, and therefore the values need to be in time since day 1 of your data.  It is recommended to uniformly distribute a value between 0 and 1 to add to your data to ensure the transformation of data is at random. 

When reading in the data, many time the first column will result in a sequence of ID numbers from 1:nrow(data), therefore the code is written to delete this first column.  If this is not the case with your data, remove this command in line 20. Should then be "data=data[idx,]".  In addition, the sample size from your data, k, can be changed to whatever is suitable.  In this demo, the full dataset only includes 200 individuals.

The remaining code specifies the length of MCMC chains, and the number of chains.  If convergence and mixing of MCMC chains is not accomplished, it is recommended to increase the length of the chains, or iteration number. 

## If recovery data is not available ...
Follow these same steps, but now use the files beginning with "oneDSA" to run the simulation. Results on recovery density will not be produced. 
