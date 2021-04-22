# Dynamic Survival Analysis Instructions
The dataset used for these demos is labled drc1_demo.csv. This is a demo dataset containing individual level infection and recovery time data can then be used in Dynamical Survival Analysis by running the following r scripts below.  The simulations are done using rstan. More information on rstan and intructions on downlaoding the package can be found here: https://mc-stan.org/users/interfaces/rstan. 

Important notes: It is highly recommended to use R Studio, as issues may occurr if not. If you run into an error ater sourcing a script, try restarting your RStudio and starting with a clean global environment.

## If recovery data is available ...
In the case that you have data for individual recovery time, as well as infection time, use the script and plot files beginning with "DSA".  Read in the data you would like to use, with the onset time first, followed by the column for recovery time. Note that the times for infection and recovery are continuous data, and therefore the values need to be in time since day 1 of your data.  It is recommended to uniformly distribute a value between 0 and 1 to add to your data to ensure the transformation of data is at random. 

When reading in the data, many times the first column will result in a sequence of ID numbers from 1:nrow(data), therefore the code is written to delete this first column.  If this is not the case with your data, remove this command in line 20. Should then be "data=data[idx,]".  In addition, the sample size from your data, k, can be changed to whatever is suitable.  In this demo, the full dataset only includes 200 individuals, so 100 of these are sampled.

The remaining code specifies the length of MCMC chains, and the number of chains.  If convergence and mixing of MCMC chains is not accomplished, it is recommended to increase the length of the chains, or iteration number. 

To run the simulation, simply source the script file, making sure the respective plot file is in the same working directory.  The resulting plots and tables will be automatically produced.

## If recovery data is not available ...
Follow these same steps, but now use the files beginning with "oneDSA" to run the simulation. Results on recovery density will not be produced. 
