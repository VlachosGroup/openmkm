This folder contains an example CTI file where lateral interactions are defined. 
After defining lateral interactions in CTI file, use the *ctml_writer.py* script 
in *<OpenMKM_Root>/scripts/* folder to convert the CTI file to XML format.

Lateral Interaction definition requires two specifications in the CTI file.
1. The surface phase has to be specified with *interacting_interface* class and
supply *interactions="all"* as one of the arguments to the class initialization.
2. Each lateral interaction defintion is specified with *lateral_interaction* keyword.
The format is *lateral_interaction("species1 species2", [Slopes list], [coverage intervals list])*.
Here *species1* is the affected species, and *species2* is the affecting species.
The *[slopes list]* contains the interaction strengths as a list. 
The *[coverage intervals list]* contains the coverage intervals for each of the 
interaction strength as a list. This requires the length of coverage intervals list to be 
length of slopes list + 1. 
However, there are some special cases. If 0 and 1 are not specified in the coverage interval 
list, they are automatically added. The previously mentioned  lengths condition is checked 
after adding 0 and 1 to the intervals list. So the definition 
*lateral_interaction('N(S1) N(S1)', [-47.0179], [0])* is valid, because the coverage 
intervals list becomes [0, 1] and its length is 2.
