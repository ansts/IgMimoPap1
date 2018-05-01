# IgMimoPap1

# Scripts, functions and data relevant to results presented on Figure 5 and Suppl. Fig. 2 - 
# comparison between IgM reactivity with different peptide libraries

# ODXsel.Rda	Big mimotope library with # copies found of each sequence  (out of 1.1x10^6 reads)
# ODXs_790_30.rda	- Mimotopes organized in 790 clusters
# aaPh7_ok.rda	- Background frequencies of aa in Ph.D-7
# mapc10post.csv	and mapc10pre.csv	- maps of the microarray gpr files
# odx79030pwmx.rda	- PWMs of 790 
# peptideLibraries.rda	
# c10.zip,	c10_class.zip and proc_c10_class.zip	- Folders with raw or processed gpr files
# zsc.csv - z1-z5 aa biophysical parameters 	

# script0.R -	Initial steps in reading the raw sequences and preparation for the Gibb's Cluster algorithm
# LOPWMtab.R, LOpwm.R, cnstrPWM.R	and cnstrPWMtab0.R -  Calculate (vectorized) PWMs of 790 clusters of mimotopes


# algrtm_c10_1.R	-  Main script (c10 pipeline) for comparison between mimotope libraries based on 10 patients' serum IgM reactivity 
# figmakerD0.R, figmakerD0_1.R	- Figures 4,6,7,8; suppl 2,4,5

#filepr.R	c10, findallscore.R,	findmaxscore.R,	findmxscore.R,	findpvscore.R,	gprLNSVM.R, gpraddNames.R, profchk.R, pscnt.R, rndpeppre.R
# - Auxiliary functions

...
