#Implementation of "Superstabilizing Control of Discrete-Time ARX Models under Error in Variables"

## Comments
in paper t = 1-na:T for y,     t = 1-nb:T-1 for u
in code t = 1:T+na for y,     t = 1:T+nb-1 for u

The design is in lambda = 1/z. Given system (Gz) in z, we first tranform it to system (Gl) in lambda, then design a data-driven superstabilizing controller (Cl) with an induced closed-loop system (ACL_l). 

The lambda quantities are transformed into the z-domain to get Cz, ACL_z

The coefficient order (na, nb) is low to high in paper (in lambda without constant term), while in MATLAB, it is high to low with constant term (in z). We performed flip, scaling and padding operations to perform the z-lambda transformations equalize them. 

## Routines
Algorithm 1
test_Full_SS.m      			 Full method for superstability, only works for low order system due to high complexity

Algorithm 2&3
test_Alt_SS.m   		         Alternatives method for superstability, works for higher order systems. Greatly reduces the complexity as compared to Full. 

test_mbc.m                               compute the model-based (ground-truth) result for comparison
complexity.m                             computes the size of psatz constraints

## Contact
Contact Tianyu (dai.ti@northeastern.edu) for any questions regarding the code
