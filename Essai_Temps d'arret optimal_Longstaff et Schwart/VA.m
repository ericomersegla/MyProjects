function Z=VA(kt,N)
 seed=rng %save seed
Z=randn(kt,N); %%v.a gaussienne 
rng(seed) %get saved seed