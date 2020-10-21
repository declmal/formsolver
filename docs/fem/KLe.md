

##### Jacobi Matrix


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
= h_{n,i} \ {}^{0}X_{in}
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
{}_{0}h_{n,j}
= \frac{\partial \ h_{n}}{\partial \ {}^{0}x_{j}}
= {}^{0}J_{jo}^{-1} \frac{\partial \ h_{n}}{\partial \ r_{o}}
= {}^{0}J_{jo}^{-1} h_{n,o}
$$

##### Displacement Increment Derivative with Respect to Global Coordinate

$$
u_{i}
= h_{n} \ U_{in}
$$

$$
{}_{0}u_{i,j}
= \frac{\partial u_{i}}{\partial \ {}^{0}x_{j}}
= {}_{0}h_{n,j} \ U_{in}
$$

##### Temporal Displacement Derivative with Respect to Global Coordinate

$$
{}^{t}u_{i}
= h_{n} \ {}^{t}U_{in}
$$

$$
{}_{0}^{t}u_{i,j}
= \frac{\partial \ {}^{t}u_{i}}{\partial \ {}^{0}x_{j}}
= {}_{0}h_{n,j} \ {}^{t}U_{in}
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
\Bigg[
{}_{0}h_{n,j} \ U_{in} +
{}_{0}h_{n,i} \ U_{jn} +
\Big(
{}_{0}^{t}u_{q,i} \ {}_{0}h_{n',j} +
{}_{0}h_{n',i} \ {}_{0}^{t}u_{q,j}
\Big) \
U_{qn'}
\Bigg]

\end{align}
$$

##### Virtual Linear Strain Tensor

$$
\begin{align}

\delta \ {}_{0}^{t}e_{ij}

&= {}_{0}e_{ij,pm} \ \delta U_{pm} \\

&= \frac{1}{2}
\Big(
{}_{0}h_{m,j} \ \delta_{ip} +
{}_{0}h_{m,i} \ \delta_{jp} +
{}_{0}^{t}u_{p,i} \ {}_{0}h_{m,j} +
{}_{0}h_{m,i} \ {}_{0}^{t}u_{p,j}
\Big) \
\delta U_{pm}

\end{align}
$$

##### Linear Strain Incremental Stiffness Matrices

[*Finite Element Procedures (2nd), P524, TABLE 6.2*]
$$
\iiint_{{}^{0}V} {}_{0}C_{ijrs} \ {}_{0}e_{rs} \ \delta \ {}_{0}e_{ij} \ \text{d} \ {}^{0}V =
\sum_{e=1}^{E} \iiint_{{}^{0}\hat V} {}_{0}C_{ijrs} \ {}_{0}e_{rs} \ \delta \ {}_{0}e_{ij} \ \text{d} \ {}^{0}\hat V
$$

Since ${}_{0}C_{ijrs} = {}_{0}C_{ijsr} = {}_{0}C_{jirs}$
$$
\begin{align}

{}_{0}C_{ijrs} \ {}_{0}^{t}e_{rs} \ \delta \ {}_{0}^{t}e_{ij}

&= \delta U_{pm} \
\Big(
{}_{0}h_{m,i} \ \delta_{jp} +
{}_{0}h_{m,i} \ {}_{0}^{t}u_{p,j}
\Big) \
{}_{0}C_{ijrs} \
\Big(
{}_{0}h_{n,r} \ \delta_{sq} +
{}_{0}h_{n,r} \ {}_{0}^{t}u_{q,s}
\Big) \
U_{qn}

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
k = 2|i-j| + \bigg\lceil \frac{i+j}{2} \bigg\rceil
$$

$$
l = 2|r-s| + \bigg\lceil \frac{r+s}{2} \bigg\rceil
$$

Thus
$$
\begin{align}

{}_{0}C_{ijrs} \ {}_{0}^{t}e_{rs} \ \delta \ {}_{0}^{t}e_{ij}

&= \delta U_{a} \ 
\Big(
{}_{0}B^{L0}_{ka} + {}_{0}^{t}B^{L1}_{ka}
\Big) \ 
{}_{0}C_{kl} \
\Big(
{}_{0}B^{L0}_{lb} + {}_{0}^{t}B^{L1}_{lb}
\Big) \
U_{b} \\

&= \delta\boldsymbol{U}^{T} \ 
\Big(
{}_{0}\boldsymbol{B}^{L0} + {}_{0}^{t}\boldsymbol{B}^{L1}
\Big)^{T} \ 
{}_{0}\boldsymbol{C} \
\Big(
{}_{0}\boldsymbol{B}^{L0} + {}_{0}^{t}\boldsymbol{B}^{L1}
\Big) \
\boldsymbol{U}

\end{align}
$$
where
$$
\begin{align}

{}_{0}B^{L0}_{ka} 

&= {}_{0}h_{m,i} \ \delta_{jp}

\end{align}
$$

$$
\begin{align}

{}_{0}^{t}B^{L1}_{ka} 

&= {}_{0}h_{m,i} \ {}_{0}^{t}u_{p,j}

\end{align}
$$

$$
{}_{0}\boldsymbol{B}^{L0}
= \left(\begin{matrix}
...& ...& ...& {}_{0}h_{m,0}& 0& 0& ...& ...& ... \\
...& ...& ...& 0& {}_{0}h_{m,1}& 0& ...& ...& ... \\
...& ...& ...& 0& 0& {}_{0}h_{m,2}& ...& ...& ... \\
...& ...& ...& {}_{0}h_{m,1}& {}_{0}h_{m,0}& 0& ...& ...& ... \\
...& ...& ...& 0& {}_{0}h_{m,2}& {}_{0}h_{m,1}& ...& ...& ... \\
...& ...& ...& {}_{0}h_{m,2}& 0& {}_{0}h_{m,0}& ...& ...& ... \\
\end{matrix}\right)
$$

$$
{}_{0}\boldsymbol{B}^{L1}
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


