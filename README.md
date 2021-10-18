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

## Reference
Tsukasa Fukunaga and Wataru Iwasaki. "Inverse Potts model improves accuracy of phylogenetic profiling." Bioinformatics Advances, vbab014.
