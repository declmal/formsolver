##### Initial Global Coordinate Interpolation

$$
{}^{0}x_{i} = {}^{0}X_{in} \ h_{n}
$$

##### Jacobi Matrix

[*Finite Element Procedures (2nd), P346*]
$$
{}^{0}\boldsymbol{J} = \left(\begin{matrix}
{}^{0}x_{0,0} & {}^{0}x_{0,1} & {}^{0}x_{0,2} \\
{}^{0}x_{1,0} & {}^{0}x_{1,1} & {}^{0}x_{1,2} \\
{}^{0}x_{2,0} & {}^{0}x_{2,1} & {}^{0}x_{2,2}
\end{matrix}\right)
$$

$$
{}^{0}J_{ij} 
= {}^{0}x_{i,j}
= \frac{\partial \ {}^{0}x_{i}}{\partial \ r_{j}}
= {}^{0}X_{in} \ h_{n,j}
$$

$$
\begin{align}

\det {}^{0}\boldsymbol{J} 

&= \epsilon_{ijk} \ {}^{0}x_{i,0} \ {}^{0}x_{j,1} \ {}^{0}x_{j,2} \\

&= ({}^{0}x_{1,0} \ {}^{0}x_{2,1} - {}^{0}x_{1,1} \ {}^{0}x_{2,0}) {}^{0}x_{0,2} +
({}^{0}x_{0,1} \ {}^{0}x_{2,0} - {}^{0}x_{0,0} \ {}^{0}x_{2,1}) {}^{0}x_{1,2} +
({}^{0}x_{0,0} \ {}^{0}x_{1,1} - {}^{0}x_{0,1} \ {}^{0}x_{1,0}) {}^{0}x_{2,2}

\end{align}
$$

$$
{}^{0}\boldsymbol{J}^{*} = 

\left(\begin{matrix}
{}^{0}x_{1,1} \ {}^{0}x_{2,2} - {}^{0}x_{1,2} \ {}^{0}x_{2,1} & 
{}^{0}x_{0,2} \ {}^{0}x_{2,1} - {}^{0}x_{0,1} \ {}^{0}x_{2,2} & 
{}^{0}x_{0,1} \ {}^{0}x_{1,2} - {}^{0}x_{0,2} \ {}^{0}x_{1,1} \\

{}^{0}x_{1,2} \ {}^{0}x_{2,0} - {}^{0}x_{1,0} \ {}^{0}x_{2,2} & 
{}^{0}x_{0,0} \ {}^{0}x_{2,2} - {}^{0}x_{0,2} \ {}^{0}x_{2,0} & 
{}^{0}x_{0,2} \ {}^{0}x_{1,0} - {}^{0}x_{0,0} \ {}^{0}x_{1,2} \\

{}^{0}x_{1,0} \ {}^{0}x_{2,1} - {}^{0}x_{1,1} \ {}^{0}x_{2,0} & 
{}^{0}x_{0,1} \ {}^{0}x_{2,0} - {}^{0}x_{0,0} \ {}^{0}x_{2,1} & 
{}^{0}x_{0,0} \ {}^{0}x_{1,1} - {}^{0}x_{0,1} \ {}^{0}x_{1,0}
\end{matrix}\right)
$$

$$
{}^{0}\boldsymbol{J}^{-1} = 
\frac{{}^{0}\boldsymbol{J^{*}}}{\det {}^{0}\boldsymbol{J}}
$$

##### Interpolation Derivative with Respect to Global Coordinate

$$
\frac{\partial \ h_{n}}{\partial \ {}^{0}\boldsymbol{x}}

= {}^{0}\boldsymbol{J}^{-1} \
\frac{\partial \ h_{n}}{\partial \ \boldsymbol{r}}
$$

$$
\frac{\partial \ h_{n}}{\partial \ {}^{0}x_{j}}
= {}^{0}J_{jo}^{-1} \frac{\partial \ h_{n}}{\partial \ r_{o}}
$$

$$
{}_{0}h_{n,j}
= h_{n,o} \ {}^{0}J_{oj}^{-1}
$$

##### Displacement Increment Derivative with Respect to Global Coordinate

$$
u_{i}
= U_{in} \ h_{n}
$$

$$
{}_{0}u_{i,j}
= \frac{\partial u_{i}}{\partial \ {}^{0}x_{j}}
= U_{in} \ {}_{0}h_{n,j}
$$

##### Temporal Displacement Derivative with Respect to Global Coordinate

$$
{}^{t}u_{i}
= {}^{t}U_{in} \ h_{n}
$$

$$
{}_{0}^{t}u_{i,j}
= \frac{\partial \ {}^{t}u_{i}}{\partial \ {}^{0}x_{j}}
= {}^{t}U_{in} \ {}_{0}h_{n,j}
$$

##### Linear Strain Tensor

$$
\begin{align}

{}_{0}^{t}e_{ij}

&= \frac{1}{2}
\Big(
{}_{0}u_{i,j} +
{}_{0}u_{j,i} +
{}_{0}^{t}u_{q,i} \ {}_{0}u_{q,j} +
{}_{0}u_{q,i} \ {}_{0}^{t}u_{q,j}
\Big) \\

&= \frac{1}{2}
U_{qn}
\Big(
{}_{0}h_{n,j} \ \delta_{iq} +
{}_{0}h_{n,i} \ \delta_{jq} +
{}_{0}^{t}u_{q,i} \ {}_{0}h_{n,j} +
{}_{0}h_{n,i} \ {}_{0}^{t}u_{q,j}
\Big) \\

\end{align}
$$

Thus
$$
{}_{0}^{t}e_{ij} = {}_{0}^{t}e_{ji}
$$

##### 2Virtual Linear Strain Tensor

$$
\delta \ {}_{0}^{t}e_{ij}

= \frac{\partial \ {}_{0}^{t}e_{ij}}{\partial \ U_{pm}} \
\delta U_{pm}

= \frac{1}{2}
\Big(
{}_{0}h_{m,j} \ \delta_{ip} +
{}_{0}h_{m,i} \ \delta_{jp} +
{}_{0}^{t}u_{p,i} \ {}_{0}h_{m,j} +
{}_{0}h_{m,i} \ {}_{0}^{t}u_{p,j}
\Big) \
\delta U_{pm}
$$

Thus
$$
\delta \ {}_{0}^{t}e_{ij} = \delta \ {}_{0}^{t}e_{ji}
$$

##### Constitutive Tensor

The following equation must be true
$$
{}_{0}^{t}C_{ijrs} = {}_{0}^{t}C_{ijsr} = {}_{0}^{t}C_{jirs}
$$
e.g. Isotropic elasticity

[*Finite Element Procedures (2nd), P194, TABLE 4.3*]
$$
{}_{0}^{t}C_{ijrs}

= \lambda \ \delta_{ij} \ \delta_{rs} + \mu \ \big(\delta_{ir} \ \delta_{js} + \delta_{is} \ \delta_{jr}\big)
$$

$$
\lambda = \frac{E \nu}{(1+\nu)(1-2\nu)}
$$

$$
\mu = \frac{E}{2(1+\nu)}
$$

##### Linear Strain Incremental Stiffness Matrix

[*Finite Element Procedures (2nd), P524, TABLE 6.2*]
$$
\iiint_{{}^{0}V} {}_{0}^{t}C_{ijrs} \ {}_{0}e_{rs} \ \delta \ {}_{0}e_{ij} \ \text{d} \ {}^{0}V

= \sum_{e=0}^{E-1} \iiint_{{}^{0}\hat V} {}_{0}^{t}C_{ijrs} \ {}_{0}e_{rs} \ \delta \ {}_{0}e_{ij} \ \text{d} \ {}^{0}\hat V
$$


$$
\begin{align}

\sum_{\{i,j\}} \sum_{\{r,s\}} {}_{0}^{t}C_{ijrs} \ {}_{0}^{t}e_{rs} \ \delta \ {}_{0}^{t}e_{ij}

=&
\sum_{\{i,j|i<j\}} \sum_{\{r,s|r<s\}} {}_{0}^{t}C_{ijrs} \ \Big(2 \ {}_{0}^{t}e_{rs}\Big) \ \Big(2 \delta \  {}_{0}^{t}e_{ij}\Big) + \\
& \sum_{\{i,j|i<j\}} \sum_{\{r,s|r=s\}} {}_{0}^{t}C_{ijrs} \ {}_{0}^{t}e_{rs} \ \Big(2 \ \delta \ {}_{0}^{t}e_{ij}\Big) + \\
& \sum_{\{i,j|i=j\}} \sum_{\{r,s|r<s\}} {}_{0}^{t}C_{ijrs} \ \Big(2 \ {}_{0}^{t}e_{rs}\Big) \ \delta \ {}_{0}^{t}e_{ij} + \\
& \sum_{\{i,j|i=j\}} \sum_{\{r,s|r=s\}} {}_{0}^{t}C_{ijrs} \ {}_{0}^{t}e_{rs} \ \delta \ {}_{0}^{t}e_{ij} \\



\end{align}
$$
Let
$$
a = 3m + p
$$

$$
b = 3n + q
$$

$$
k = 2(j-i) + \bigg\lceil \frac{i+j}{2} \bigg\rceil
$$

$$
l = 2(s-r) + \bigg\lceil \frac{r+s}{2} \bigg\rceil
$$

###### Routine Component of Linear Strain Displacement Matrix

[*Finite Element Procedures (2nd), P555, TABLE 6.6*]
$$
{}_{0}\boldsymbol{B}_{L0}
= \left(\begin{matrix}
...& ...& ...& {}_{0}h_{m,0}& 0& 0& ...& ...& ... \\
...& ...& ...& 0& {}_{0}h_{m,1}& 0& ...& ...& ... \\
...& ...& ...& 0& 0& {}_{0}h_{m,2}& ...& ...& ... \\
...& ...& ...& {}_{0}h_{m,1}& {}_{0}h_{m,0}& 0& ...& ...& ... \\
...& ...& ...& 0& {}_{0}h_{m,2}& {}_{0}h_{m,1}& ...& ...& ... \\
...& ...& ...& {}_{0}h_{m,2}& 0& {}_{0}h_{m,0}& ...& ...& ...
\end{matrix}\right)
$$

###### Bbar Approach

$$
{}_{0}\boldsymbol{B}_{L0}^{dil}
= \frac{1}{3} \left(\begin{matrix}
... & ... & ... & {}_{0}h_{m,0} & {}_{0}h_{m,1} & {}_{0}h_{m,2} & ... & ... & ... \\
... & ... & ... & {}_{0}h_{m,0} & {}_{0}h_{m,1} & {}_{0}h_{m,2} & ... & ... & ... \\
... & ... & ... & {}_{0}h_{m,0} & {}_{0}h_{m,1} & {}_{0}h_{m,2} & ... & ... & ... \\
... & ... & ... & 0 & 0 & 0 & ... & ... & ... \\
... & ... & ... & 0 & 0 & 0 & ... & ... & ... \\
... & ... & ... & 0 & 0 & 0 & ... & ... & ...
\end{matrix}\right)
$$

$$
{}^{0}V

= \iiint_{{}^{0}\hat V} \text{d} \ {}^{0} \hat V

= \sum_{i=0}^{\bar N_{int}-1}
\omega_{i} \
\Big|\det {}^{0} \boldsymbol{J}\Big|
$$

$$
{}_{0}\boldsymbol{\bar B}_{L0}^{dil}

= \frac{1}{{}^{0}V} \ 
\iiint_{{}^{0}\hat V} \
{}_{0}\boldsymbol{B}_{L0}^{dil} \ 
\text{d} \ {}^{0}V

= \frac{1}{{}^{0}V}
\sum_{i=0}^{\bar N_{int}-1} \
\omega_{i} \
{}_{0}\boldsymbol{B}_{L0}^{dil} \
\Big|\det {}^{0}\boldsymbol{J} \Big|
$$

$$
{}_{0}\boldsymbol{\bar B}_{L0}

= {}_{0}\boldsymbol{B}_{L0} -  {}_{0}\boldsymbol{B}_{L0}^{dil} + {}_{0}\boldsymbol{\bar B}_{L0}^{dil}
$$

###### Initial Displacement Component of Linear Strain Displacement Matrix

$$
{}_{0}^{t}\boldsymbol{B}_{L1}
= \left(\begin{matrix}
...& ...& ...& 
{}_{0}h_{m,0} \ {}_{0}^{t}u_{0,0}& 
{}_{0}h_{m,0} \ {}_{0}^{t}u_{1,0}& 
{}_{0}h_{m,0} \ {}_{0}^{t}u_{2,0}& ...& ...& ... \\

...& ...& ...& 
{}_{0}h_{m,1} \ {}_{0}^{t}u_{0,1}& 
{}_{0}h_{m,1} \ {}_{0}^{t}u_{1,1}& 
{}_{0}h_{m,1} \ {}_{0}^{t}u_{2,1}& ...& ...& ... \\

...& ...& ...& 
{}_{0}h_{m,2} \ {}_{0}^{t}u_{0,2}& 
{}_{0}h_{m,2} \ {}_{0}^{t}u_{1,2}& 
{}_{0}h_{m,2} \ {}_{0}^{t}u_{2,2}& ...& ...& ... \\

...& ...& ...& 
{}_{0}h_{m,0} \ {}_{0}^{t}u_{0,1} + {}_{0}h_{m,1} \ {}_{0}^{t}u_{0,0}& 
{}_{0}h_{m,0} \ {}_{0}^{t}u_{1,1} + {}_{0}h_{m,1} \ {}_{0}^{t}u_{1,0}& 
{}_{0}h_{m,0} \ {}_{0}^{t}u_{2,1} + {}_{0}h_{m,1} \ {}_{0}^{t}u_{2,0}& ...& ...& ... \\

...& ...& ...& 
{}_{0}h_{m,1} \ {}_{0}^{t}u_{0,2} + {}_{0}h_{m,2} \ {}_{0}^{t}u_{0,1}& 
{}_{0}h_{m,1} \ {}_{0}^{t}u_{1,2} + {}_{0}h_{m,2} \ {}_{0}^{t}u_{1,1}& 
{}_{0}h_{m,1} \ {}_{0}^{t}u_{2,2} + {}_{0}h_{m,2} \ {}_{0}^{t}u_{2,1}& ...& ...& ... \\

...& ...& ...& 
{}_{0}h_{m,0} \ {}_{0}^{t}u_{0,2} + {}_{0}h_{m,2} \ {}_{0}^{t}u_{0,0}& 
{}_{0}h_{m,0} \ {}_{0}^{t}u_{1,2} + {}_{0}h_{m,2} \ {}_{0}^{t}u_{1,0}& 
{}_{0}h_{m,0} \ {}_{0}^{t}u_{2,2} + {}_{0}h_{m,2} \ {}_{0}^{t}u_{2,0}& ...& ...& ... \\
\end{matrix}\right)
$$

###### Linear Strain Displacement Matrix

$$
{}_{0}^{t}\boldsymbol{B}_{L}
= {}_{0}\boldsymbol{\bar B}_{L0} + {}_{0}^{t}\boldsymbol{B}_{L1}
$$

##### Linear Component of Element Stiffness Matrix

$$
\boldsymbol{K}_{L}^{e}

= \iiint_{{}^{0}\hat V}
{}_{0}^{t}\boldsymbol{B}_{L}^{T} \ 
{}_{0}^{t}\boldsymbol{C} \ 
{}_{0}^{t}\boldsymbol{B}_{L} \ \text{d}{}^{0}\hat V

= \sum_{i=0}^{\hat N_{int}-1}
\omega_{i} \ 
{}_{0}^{t}\boldsymbol{B}_{L}(\boldsymbol{r}_{i}) \ 
{}_{0}^{t}\boldsymbol{C} \
{}_{0}^{t}\boldsymbol{B}_{L}(\boldsymbol{r}_{i}) \
\Big|\det {}^{0}\boldsymbol{J}\Big|
$$

##### Nonlinear Strain Tensor

$$
{}_{0}\eta_{ij}

= \frac{1}{2} {}_{0}u_{q,i} \ {}_{0}u_{q,j}

= \frac{1}{2} {}_{0}h_{n',i} \ {}_{0}h_{n,j} \ U_{qn'} \ U_{qn}
$$

##### Virtual Nonlinear Strain Tensor

$$
\begin{align}

\delta \ {}_{0}\eta_{ij}

&= {}_{0}\eta_{ij,pm} \ \delta U_{pm} \\

&= \frac{1}{2} {}_{0}h_{n',i} \ {}_{0}h_{n,j}
\Big(
U_{qn} \ \delta U_{qn'} +
U_{qn'} \ \delta U_{qn}
\Big)

\end{align}
$$

Thus
$$
\delta \ {}_{0}\eta_{ij} = \delta \ {}_{0}\eta_{ji}
$$

#### Second Piola-Kirchoff Stress Tensor

The following equation must be true
$$
{}_{0}^{t} S_{ij} = {}_{0}^{t} S_{ji}
$$
e.g. 
$$
{}_{0}^{t}S_{ij} = {}_{0}^{t}C_{ijrs} \ {}_{0}^{t}\epsilon_{rs}
$$

##### Nonlinear Strain Incremental Stiffness Matrix

[*Finite Element Procedures (2nd), P524, TABLE 6.2*]
$$
\iiint_{{}^{0}V} {}_{0}^{t}S_{ij} \ \delta \ {}_{0}\eta_{ij} \ \text{d} \ {}^{0}V =

\sum_{e=0}^{E-1} \iiint_{{}^{0}\hat V} {}_{0}^{t}S_{ij} \ \delta \ {}_{0}\eta_{ij} \ \text{d} \ {}^{0}\hat V
$$

$$
\begin{align}

{}_{0}^{t}S_{ij} \ \delta \ {}_{0}\eta_{ij}

&= \delta U_{pm} \ {}_{0}h_{m,i} \ \delta_{pr} \ {}_{0}^{t}S_{ij} \ \delta_{rs} \ {}_{0}h_{n,j} \ \delta_{qs} \ U_{qn}

\end{align}
$$
Let
$$
a = 3m + p
$$

$$
b = 3n + q
$$

$$
k = 3r + i
$$

$$
l = 3s + j
$$

Thus
$$
\begin{align}

{}_{0}^{t}S_{ij} \ \delta \ {}_{0}\eta_{ij}

&= \delta U_{a} \ 
{}_{0}B^{NL}_{ka} \ {}_{0}^{t}S_{kl} \ {}_{0}B^{NL}_{lb} \ U_{b} \\

&= \delta\boldsymbol{U}^{T} \ 
\Big(
{}_{0}\boldsymbol{B}^{NL}
\Big)^{T} \ 
{}_{0}^{t}\boldsymbol{S} \ {}_{0}\boldsymbol{B}^{NL} \ \boldsymbol{U}

\end{align}
$$

$$
{}_{0}\boldsymbol{C} = \left(\begin{matrix}
\lambda + 2\mu & \lambda & \lambda & 0 & 0 & 0 \\
\lambda & \lambda + 2\mu & \lambda & 0 & 0 & 0 \\
\lambda & \lambda & \lambda + 2\mu & 0 & 0 & 0 \\
0 & 0 & 0 & \mu & 0 & 0 \\
0 & 0 & 0 & 0 & \mu & 0 \\
0 & 0 & 0 & 0 & 0 & \mu
\end{matrix}\right)
$$



#### Incompatible

$$
u_{i}
= U_{in} \ h_{n} + 
U_{in}^{*} \ h_{n}^{*}
$$

$$
{}_{0}u_{i,j}
= \frac{\partial u_{i}}{\partial \ {}^{0}x_{j}}
= U_{in} \ {}_{0}h_{n,j} + U_{in}^{*} \ {}_{0}h_{n,j}^{*}
$$

$$
\begin{align}

{}_{0}e_{ij}

&= \frac{1}{2}
\Big(
{}_{0}u_{i,j} +
{}_{0}u_{j,i}
\Big) \\

&= \frac{1}{2}
U_{qn}
\Big(
{}_{0}h_{n,j} \ \delta_{iq} +
{}_{0}h_{n,i} \ \delta_{jq}
\Big) +
\frac{1}{2}
U_{qn}^{*}
\Big(
{}_{0}h_{n,j}^{*} \ \delta_{iq} +
{}_{0}h_{n,i}^{*} \ \delta_{jq}
\Big)\\

\end{align}
$$

$$
{}_{0}e_{ij} = {}_{0}e_{ji}
$$




$$
\begin{align}

\delta \ {}_{0}e_{ij}

&= \frac{\partial \ {}_{0}e_{ij}}{\partial \ U_{pm}} \
\delta U_{pm} + 
\frac{\partial \ {}_{0}e_{ij}}{\partial \ U_{pm}^{*}} \
\delta U_{pm}^{*} \\

&= \frac{1}{2}
\Big(
{}_{0}h_{m,j} \ \delta_{ip} +
{}_{0}h_{m,i} \ \delta_{jp}
\Big) \
\delta U_{pm} +
\frac{1}{2}
\Big(
{}_{0}h_{m,j}^{*} \ \delta_{ip} +
{}_{0}h_{m,i}^{*} \ \delta_{jp}
\Big) \
\delta U_{pm}^{*}
\end{align}
$$

$$
\delta \ {}_{0}e_{ij} = \delta \ {}_{0}e_{ji}
$$


$$
{}_{0}^{t}C_{ijrs} \ {}_{0}e_{rs} \ \delta \ {}_{0}e_{ij}

= 
$$
