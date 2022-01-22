# Ipm
Ipm is a program for calculating direct information based on the inverse Potts Model using the persistative contrastive divergence method.

## Version
Version 0.1.0 (2021/10/18)

## Options
    (Required)
        -i STR    InputFileName
        -o STR    OutputFileName
        
    (Optional) 
        -e DBL    The learning rate in parameter updates. [default:0.01]
        -l DBL    The hyperparameter of l2-regularization. [default: 0]
        -c INT    The number of parameter updates in Gibbs sampling. [default: 3000]
        -s INT    The sample size in Gibbs sampling. [default: 200]

## Input file format
"n l q" is described in the first line. The second line specifies the number of the genome, the number of gene families, and the number of states. After the third line, the state of each gene family in each genome is described.

## Example of input file format
    n l q
    3 5 2
    1 0 0 1 1
    1 0 1 0 1
    0 1 0 1 0


## Reference
Tsukasa Fukunaga and Wataru Iwasaki. "Inverse Potts model improves accuracy of phylogenetic profiling." in press, Bioinformatics
