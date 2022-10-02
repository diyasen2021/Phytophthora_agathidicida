This example will create sum of values over bins in R using dplyr.

Read in data file, repeat_withlength_oct2022.txt. This is output from repeatmodeller
and repeatmasker. Length is number of bases of each repeat, calculated by End - Start columns. 

To sum over length for bins of specific length, binlength.R uses common dplyr functions such as 
group_by, mutate and distinct.  
