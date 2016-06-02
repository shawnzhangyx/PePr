
> PePr v1.1.10 
> Timestamp: 06-01-2016 Yanxiao Zhang zhangyx.shawn@gmail.com
------------------------------------------------------------

###Introduction
PePr is a ChIP-Seq Peak-calling and Prioritization pipeline 
that uses a sliding window approach and models read counts across 
replicates and between groups with a negative binomial distribution. 
PePr empirically estimates the optimal shift/fragment size and 
sliding window width, and estimates dispersion from the local genomic
area. Regions with less variability across replicates are ranked more
favorably than regions with greater variability. Optional 
post-processing steps are also made available to filter out peaks
not exhibiting the expected shift size and/or to narrow the width of peaks.

### Installation
1. Make sure your python version is higher than 2.6. Version 3.X may not be fully supported yet.
2. Install **pip** in your system if you don't have it. 
3. `pip install PePr` or `pip install PePr --user`(if you don't have administrator privilege). Optionally, you can download tarball (PePr-[version].tar.gz) from [github](https://github.com/zhangyx.shawn/PePr/) and install using `pip install PePr-[version].tar.gz`
4. If installation is successful, you could directly invoke the script by typing `PePr`. A help message will show up. 

### Supported File Formats
* Single-end: BED, BAM, SAM. 
* Paired-end: BAM, SAM. The files must be sorted by the read names. Users can use `samtools sort -n sample.bam sample.sorted_by_name` to sort the file. 

### Basic Usage Examples
*Warning: These are working examples with minimal required parameters. For the best performance (or to avoid bad fitting) on your data, please read the manual carefully and choose the right parameters.* 
* For peak-calling, run: `PePr -c chip_rep1.bam,chip_rep2.bam -i input_rep1.bam,input_rep2.bam -f bam`
* For differential binding analysis with input samples, run: `PePr -c chip1_rep1.bam,chip1_rep2.bam -i input1_rep1.bam,input1_rep2.bam --chip2 chip2_rep1.bam,chip2_rep2.bam --input2 input2_rep2.bam,input2_rep2.bam -f bam`
* For differential binding analysis without input samples, run: `PePr -c chip1_rep1.bam,chip1_rep2.bam --chip2 chip2_rep1.bam,chip2_rep2.bam -f bam`
  

###LINKS
* https://github.com/zhangyx.shawn/PePr/ # Source code
* https://ones.ccmb.med.umich.edu/tags/PePr/ # PePr FAQ
* https://pypi.python.org/pypi/pepr #PyPI package index



### QUESTIONS?
Questions are preferred to be posted on https://ones.ccmb.med.umich.edu/, you'll very likely find similar problems in there too (https://ones.ccmb.med.umich.edu/tags/PePr/?tab=Summary). You're also welcome to shoot me an e-mail at yanxiazh@umich.edu, I'll try replying to you as soon as possible. 


###CITATION
Zhang Y, Lin YH, Johnson TD, Rozek LS, Sartor MA. [PePr: A peak-calling prioritization pipeline to identify consistent or differential peaks from replicated ChIP-Seq data.](http://www.ncbi.nlm.nih.gov/pubmed/24894502) Bioinformatics. 2014.
