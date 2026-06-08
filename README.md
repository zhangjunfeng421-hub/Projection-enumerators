### Projection-enumerators
This project includes all R codes and results about the simulation and examples of a statistical article, "A General Theory of Projection Enumerators for Factorial and Space-Filling Designs". We illustrate each file next.

## PSE.R
This code mainly includes four useful functions for a given design D with levels $s^p$. One is "PSE(D,s,p,x,y)", which is used to compute $E_{\chi}(D;x,y)$. To obtain the corresponding values of the truncated enumerator $E^(k)_{\chi}(D;x,y)$, one should apply $TSPE(D,s,p,x,y,t)$.
PSFPI(D,s,p) is used to obtain all $P_{i,j}(D)$ fast, and TPSFPI(D,s,p,t) is used to obtain $P_{i,j}(D)$ ($j \leq t$) fast.
