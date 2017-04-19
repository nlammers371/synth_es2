-------------OVERVIEW----------------

The goal of this project is to explore the importance of Eve Stripe Two sequence organization at multiple scales. 

Two alternative metrics were explored as potential metrics for the importance of regions within the enhancer to its overall function:

--Conservation across related species
--Presence of Predicted/consensus TFBS

For each metrics, "meaningful" regions were identified and the full enhancer sequence was segmented into subsections that preserved the contiguity of these regions.
Region position was then shufled such that no segment remained in original position and no segment had original neighbors. Thus both relative and global
structure were perturbed. A hybrid algorithmic/stochastic approach was used to search for novel arrangements that minimized disruption of segment boundaries.
Shuffling procedure was employed for regions of 6 or greater and 32 or greater bp for both conserved and TFBS schemes to probe relevance of scale.

--------------Code------------------ 
.code\analysis\
Files used for core of analysis are in .code\analysis\ directory. Numbers indicate order in which files should be run. Files without number are either exploratory in nature or checks.

.code\functions\
contains all relevant function files used in analysis

.code\imports\
relevant import files.

------------Data-----------------
Data files used for analysis can be found in \data\analysis\ToFIMO\Results_12.01.16

