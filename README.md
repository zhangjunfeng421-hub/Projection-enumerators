### Projection-enumerators
This project includes all R codes and results about the simulation and examples of a statistical article, "A General Theory of Projection Enumerators for Factorial and Space-Filling Designs". We illustrate each file next.

## PSE.R
This code mainly includes four useful functions for a given design $D$ with levels $s^p$. One is "PSE(x,y,D,s,p)", which is used to compute $E_{\chi}(D;x,y)$. To obtain the corresponding values of the truncated enumerator $E^{(k)}_{\chi}(D;x,y)$, one should apply "TSPE(x,y, D,s,p,t)".
"PSFPI(D,s,p)" is used to obtain all $P_{i,j}(D)$ fast, and "TPSFPI(D,s,p,t)" is used to obtain $P_{i,j}(D)$ ($j \leq t$) fast.

## PBE.R
This parallels "PSE.R", and it is designed for $B_{i,j}(D)$, where the design $D$ has $s$ levels. "PBE(x,y,D,s)" is used to compute $E_{\varphi}(D;x,y)$. To obtain the corresponding values of the truncated enumerator $E^(k)_{\varphi}(D;x,y)$, one should apply "TPBE(x,y, D,s,t)".
"PWPI(D,s)" is used to obtain all $B_{i,j}(D)$ fast, and "TPWPI(D,s,t)" is used to obtain $B_{i,j}(D)$ ($j \leq t$) fast.

## GSOA_opt.R
This code file only has one useful function, "GSOA_opt(polyprimitive, s, K = (s - 1), q = (length(polyprimitive)-1), add = 0)" to construct the optimal GSOA(s^p,K*(s^p-1)/(s-1),s^q,t). We need to illustrate the inputs "polyprimitive", "q",  and "add". 

"polyprimitive" is a coefficient vector to determine the primitive polynomial. For instance, if we are going to construct a design with $3^2$ levels with the primitive polynomial $x^2+x+2$, then the vector should be $c(1,1,2)$. If we are going to construct a design with $3^3$ levels with the primitive polynomial $x^3+2x+1$, then the vector should be $c(1,0,2,1)$.

$q$ is a number input, so the design levels are "0,\ldots,s^q-1". In the default setting, the number of levels is $s^p$. For example, with the primitive polynomial $x^3+2x+1$ and $s=3$, the number of levels should be $3^3=27$; but setting $q=2$ generates an optimal GSOA with $3^2=9$ levels.

Tian and Xu (2024) claim that the shift of the 
