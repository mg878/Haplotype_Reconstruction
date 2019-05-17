# A Novel Haplotype Reconstruction Method (MLHapRec)
`MLHapRec` is a maximum likelihood approach for haplotype reconstruction which uses the [SAMFIRE](https://github.com/cjri/samfire/blob/master/README.md) multi-locus variant calling package for processing short-read data where samples are collected during a transmission of a virus (like influenza) from a donor and a recipient population (i.e. two time points). 

This code is originally developed for estimating the transmission bottleneck size of flu, but can be used for any other set of short-read data for the purpose of constructing the underlying (unkonw) haplotypes and estimating bottleneck size. 

## Requirements
This code requires a `Multi_locus_trajectories.out` file and an inferred noise parameter `C` from `SAMFIRE`. An initial seed value and number of attempts for reconstruction can also be provided, otherwise a seed = 1 and attempts = 20 are used. Number of attempts determine, for every time a new haplotype is added to the list, how many times a set of candidate haplotypes are optimised before the code accepts the list and compares the BIC value with the previous step. 

## Usage
Given that you have inferred the noise parameter, e.g. C = 660, you can run the code by typing:
> ./run_MLHapRec /path/to/directory/Multi_locus_trajectories.out 660 123 10

where we used a seed value 123 and 10 attempts.

## Output
Running the code, typically, takes anywhere from 10sec to 10h depending on the size of `Multi_locus_trajectories.out` and number of attempts taken. For instance, if there are about 10 partial haplotype reads in the `Multi_locus_trajectories.out` and 20 attempts are made for reconstruction each time, it takes up to ~15min for the full haplotype reconstruction. 

After a successful execution, `MLHapRec` generates the following **three** files: 
- `update.txt` shows how the optimisation process takes place at every round as a new haplotype is added to the list of inferred haplotypes 
- `raw_haplotypes.txt` lists all the reconstructed haplotypes with their frequencies before and after transmission.
- `outcome_1.txt` lists the haplotypes from `raw_haplotypes.txt` which have not acquired a mutation once established in the recipient -- a mutation is defined to be a haplotype with frequency <1e-7 in the donor and >1e-7 in the recipient.
Note that in the case of no mutation, `raw_haplotypes.txt` and `outcome_1.txt` are identical.
