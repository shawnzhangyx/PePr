
> Timestamp: 03-15-2015 Yanxiao Zhang zhangyx.shawn@gmail.com
------------------------------------------------------------

###INTRODUCTION
PePr is a ChIP-Seq Peak-calling and Prioritization pipeline 
that uses a sliding window approach and models read counts across 
replicates and between groups with a negative binomial distribution. 
PePr empirically estimates the optimal shift/fragment size and 
sliding window width, and estimates dispersion from the local genomic
area. Regions with less variability across replicates are ranked more
favorably than regions with greater variability. Optional 
post-processing steps are also made available to filter out peaks
not exhibiting the expected shift size and/or to narrow the width of peaks.

###LINKS
* https://github.com/troublezhang/PePr/ # Source code
* https://ones.ccmb.med.umich.edu/wiki/PePr/ # Manual
* https://pypi.python.org/pypi/pepr #PyPI package index


###INSTALL & USAGE
Please visit our wiki page(https://ones.ccmb.med.umich.edu/wiki/PePr/) for instructions. 

### QUESTIONS?
Questions are preferred to be posted on https://ones.ccmb.med.umich.edu/, you'll very likely find similar problems in there too (https://ones.ccmb.med.umich.edu/tags/PePr/?tab=Summary). You're also welcome to shoot me an e-mail at yanxiazh@umich.edu, I'll try replying to you as soon as possible. 


###CITATION
Zhang Y, Lin YH, Johnson TD, Rozek LS, Sartor MA. PePr: A peak-calling prioritization pipeline to identify consistent or differential peaks from replicated ChIP-Seq data. Bioinformatics. 2014.
