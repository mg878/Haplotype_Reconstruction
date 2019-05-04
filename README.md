# A Novel Haplotype Reconstruction Method (MLHapRec)
`MLHapRec` is a maximum likelihood approach for haplotype reconstruction which uses the [SAMFIRE](https://github.com/cjri/samfire/blob/master/README.md) multi-locus variant calling technique for processing short-read data where samples are collected during a transmission of a virus (like influenza) from a donor and a recipient population (i.e. two time points). 

This code is originally developed for estimating the transmission bottleneck size of flu, but can be used for any other set of short-read data for the purpose of constructing the underlying (unkonw) haplotypes. 

## Requirements
This code requires `Multi_locus_trajectories.out` file and an inferred noise parameter `C` from `SAMFIRE`.

## Usage
Given that you have inferred the noise parameter, e.g. C = 660, you can run the code by typing:
> ./run_MLHapRec /path/to/directory/Multi_locus_trajectories.out 660

## Output
Running the code, typically, takes anywhere from 10sec to 10h depending on the size of `Multi_locus_trajectories.out` -- If there are about 10 partial haplotype reads in `Multi_locus_trajectories.out`, it takes up to ~15min for the full haplotype reconstruction. Note that the nature of haplotype reconstruction using `MLHapRec` is stochastic. Therefore, not every run would give **identical** set of reconstructed haplotypes. 

After a successful execution, `MLHapRec` generates the following **three** files: 
- `update.txt` shows how the optimisation process takes place at every round as a new haplotype is added to the list of inferred haplotypes 
- `raw_haplotypes.txt` lists all the reconstructed haplotypes with their frequencies before and after transmission.
- `outcome_1.txt` lists the haplotypes from `raw_haplotypes.txt` which have not acquired a mutation once established in the recipient.

## Description
Here, we discuss how `MLHapRec` process the data and makes the haplotype set:

