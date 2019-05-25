Chemkin files are widely used to distribute reaction mechanisms. Also some times, 
they are easier to work with.  [Cantera website](https://cantera.org/tutorials/input-files.html) 
has excellent guide on converting chemkin files to CTI and ultimately to XML format. 
However the script used for conversion can fail some times. Here are some steps you 
can take to overcome those issues.

1. Try to use  *<OpenMKM_Root>/scripts/ck2htrct.py* in place of ck2cti .
2. Read the error messages and remove the failing parts  and rerun. Add the missing parts manually. 
3. Lateral interactions can be input eiter in CTI/XML file. Use CTI to specify the lateral interactions
and use the supplied script *<OpenMKM_Root>/scripts/ctmlwriter.py* to get the XML file.

