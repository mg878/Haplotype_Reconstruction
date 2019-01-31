# A Novel Haplotype Reconstruction Method (MLHapRec)
MLHapRec is a maximum likelihood approach for haplotype reconstruction which uses the [SAMFIER](https://github.com/cjri/samfire/blob/master/README.md) multi-locus variant calling technique for processing short-read data where samples are collected during a transmission of a virus (like influenza) from a donor to a recipient population. 

This code is originally developed for estimating the transmission bottleneck size of flu, but can be used for any other set of short-read data. Below, we will discuss a scenario where we have 100 transmission events for flu virus collected from two time-points (before and after transmission). 

## Format of the data
This code requires a file(s) which contains all the 8 segments of flu in the following order: HA, MP, NA, NP, NS, PA, PB1, and PB2 
Within each segment, the code looks for the following output files from SAMFIER: Multi_locus_trajectories.out, Loci*.dat, and Hap_data*.dat


## Compiling the code

1- **open** blah
- `print("Hello")` 
- hi  
