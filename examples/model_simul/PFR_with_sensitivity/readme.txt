Plug flow reactor is modeled and sensitivity coefficients are calculated simultaneously.
 
There are two approaches to sensitivity.

1) Fisher Information Matrix (FIM)
Here the objective is to prune out unimportant reactions from a large set of reactions. Refer to the input file `pfr_1d.yaml` in FIM folder on how to enable FIM. The computed diagonal elements of Fisher Information Matrix are printed in 1d_pfr_sensitivity file. Reactions with high FIM values may be important or not, but those with low values are definitely not. Using the low FIM values as a guide, one can prune out the unimportant values.

2) Local sensitivity analysis (LSA)
This step can be performed after doing FIM or can be done as stand alone step.
The reactions for which sensitivity coefficients are calculated are listed under `sensitivity` keyword.

The output file 1d_pfr_sensitivity.out lists the sensitivity coefficients of the state variables, gas mass fractions and surface coverages w.r.t. the reactions. Each row corresponds to one reaction and the columns represent the state variables, ...

If the values are separated by commas, use pandas.read_csv function to load the data into python.
