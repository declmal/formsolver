# Doctoral Thesis Framework

[TOC]

$$
\boldsymbol{u}_{\text{xfem}}(\boldsymbol{x}) =
\sum_{j=1}^{n} N_{j}(\boldsymbol{x}) \boldsymbol{u}_{j}
+ \sum_{h=1}^{mh} N_{h} \Big( H(x) - H_{h}(x_{h}) \Big) \boldsymbol{a}_{h}
+ \sum_{k=1}^{mt} N_{k}(\boldsymbol{x}) \Bigg[ \sum_{l=1}^{4} \Big( F_{l}(\boldsymbol{x}) - F_{l}(\boldsymbol{x}_{k}) \Big) \boldsymbol{b}_{k}^{l} \Bigg]
$$

## CHAPTER 2

### Iterator

load increment $\Delta t$

Convergence Criteria

Displacement Improver

#### Full Newton-Raphson Iterator

1. Beginning of Increment

Initial conditions. [*Finite Element Procedures (2nd), P755*]
$$
{}^{t + \Delta t}\boldsymbol{K}^{(0)} = {}^{t}\boldsymbol{K}
$$

$$
{}^{t + \Delta t}\boldsymbol{F}^{(0)} = {}^{t}\boldsymbol{F}
$$

$$
{}^{t + \Delta t}\boldsymbol{U}^{(0)} = {}^{t}\boldsymbol{U}
$$

Evaluate ${}^{t + \Delta t}\boldsymbol{R}$.

2. Each Iteration

Evaluate nodal point force vectors ${}^{t + \Delta t}\boldsymbol{F}^{(i-1)}$.
$$
{}^{t + \Delta t}\boldsymbol{F}^{(i-1)}
= \boldsymbol{F} \Big|_{{}^{t + \Delta t} \boldsymbol{U}^{(i-1)}}
$$
Evaluate tangent stiffness matrix ${}^{t + \Delta t}\boldsymbol{K}^{(i-1)} $  [*Finite Element Procedures (2nd), P755*], it should be consistent with ${}^{t + \Delta t}\boldsymbol{F}^{(i-1)}$ [*Finite Element Procedures (2nd), P757-758*] and nonsingular [*Finite Element Procedures (2nd), P757*] for convergence.
$$
{}^{t + \Delta t}\boldsymbol{K}^{(i-1)}
= \frac{\partial \boldsymbol{F}}{\partial \boldsymbol{U}} \Bigg|_{{}^{t + \Delta t} \boldsymbol{U}^{(i-1)}}
$$
Solve for $\Delta \boldsymbol{U}^{(i)}$. [*Finite Element Procedures (2nd), P755, (8.88)*]
$$
{}^{t + \Delta t}\boldsymbol{K}^{(i-1)} \Delta \boldsymbol{U}^{(i)} 
= {}^{t + \Delta t}\boldsymbol{R} - {}^{t + \Delta t}\boldsymbol{F}^{(i-1)}
$$

Search for the improved solution ${}^{t + \Delta t}\boldsymbol{U}^{(i)}$.

3. End of Increment

Store:
$$
{}^{t + \Delta t}\boldsymbol{K}
$$

$$
{}^{t + \Delta t}\boldsymbol{F}
$$

$$
{}^{t + \Delta t}\boldsymbol{U}
$$

#### Modified Newton-Raphson Iterator

1. Beginning of Increment

Initial conditions [*Finite Element Procedures (2nd), P759*].
$$
{}^{t + \Delta t}\boldsymbol{F}^{(0)} = {}^{t}\boldsymbol{F}
$$

$$
{}^{t + \Delta t}\boldsymbol{U}^{(0)} = {}^{t}\boldsymbol{U}
$$

Determine whether to update tangent stiffness matrix [*Finite Element Procedures (2nd), P759*] or use other self-adaptive ways.
$$
{}^{\tau}\boldsymbol{K} =
\frac{\partial \boldsymbol{F}}{\partial \boldsymbol{U}} \Bigg|_{{}^{t} \boldsymbol{U}}
$$
Evaluate ${}^{t + \Delta t}\boldsymbol{R}$.

2. Each Iteration

Evaluate nodal point force vectors ${}^{t + \Delta t}\boldsymbol{F}^{(i-1)}$.
$$
{}^{t + \Delta t}\boldsymbol{F}^{(i-1)}
= \boldsymbol{F} \Big|_{{}^{t + \Delta t} \boldsymbol{U}^{(i-1)}}
$$
Solve for $\Delta \boldsymbol{U}^{(i)}$. [*Finite Element Procedures (2nd), P759, (8.93)*]
$$
{}^{\tau}\boldsymbol{K} \Delta\boldsymbol{U}^{(i)} 
= {}^{t + \Delta t}\boldsymbol{R} - {}^{t + \Delta t}\boldsymbol{F}^{(i-1)}
$$
Search for the improved solution ${}^{t + \Delta t}\boldsymbol{U}^{(i)}$.

3. End of Increment

Store:
$$
{}^{t + \Delta t}\boldsymbol{F}
$$

$$
{}^{t + \Delta t}\boldsymbol{U}
$$

$$
{}^{\tau}\boldsymbol{K}
$$

#### BFGS Iterator

1. Beginning of Increment

Determine whether to update tangent stiffness matrix or use other self-adaptive ways.
$$
{}^{\tau}\boldsymbol{K}^{-1} =
\Bigg(\frac{\partial \boldsymbol{F}}{\partial \boldsymbol{U}} \Bigg|_{{}^{t} \boldsymbol{U}}\Bigg)^{-1}
$$
Evaluate ${}^{t + \Delta t}\boldsymbol{R}$.

Initial conditions [*Finite Element Procedures (2nd), P754, (8.80)*] and [*Finite Element Procedures (2nd), P760*]
$$
\Delta\boldsymbol{R}^{(0)} = {}^{t + \Delta t}\boldsymbol{R} - {}^{t}\boldsymbol{F}
$$

$$
\Big({}^{t + \Delta t}\boldsymbol{K}^{-1}\Big)^{(0)} = {}^{\tau}\boldsymbol{K}^{-1}
$$

2. Each Iteration

Displacement vector Increment Evaluation [*Finite Element Procedures (2nd), P760, (8.97)*].
$$
\Delta \boldsymbol{U}^{(i)} =
\Big({}^{t + \Delta t}\boldsymbol{K}^{-1}\Big)^{(i-1)}
\Delta\boldsymbol{R}^{(i-1)}
$$
Search for the improved solution ${}^{t + \Delta t}\boldsymbol{U}^{(i)}$.

Displacement Increment Calculation [*Finite Element Procedures (2nd), P759, (8.94)*].
$$
\boldsymbol{\delta}^{(i)} = {}^{t + \Delta t}\boldsymbol{U}^{(i)} - {}^{t + \Delta t}\boldsymbol{U}^{(i-1)}
$$
Evaluate nodal point force vectors ${}^{t + \Delta t}\boldsymbol{F}^{(i)}$.
$$
{}^{t + \Delta t}\boldsymbol{F}^{(i)}
= \boldsymbol{F} \Big|_{{}^{t + \Delta t} \boldsymbol{U}^{(i)}}
$$
Out of Balance Load Calculation [*Finite Element Procedures (2nd), P754, (8.80)*].
$$
\Delta\boldsymbol{R}^{(i)} = {}^{t + \Delta t}\boldsymbol{R} - {}^{t + \Delta t}\boldsymbol{F}^{(i)}
$$
Increment of Out-of-balance Loads Calculation [*Finite Element Procedures (2nd), P759, (8.95)*].
$$
\boldsymbol{\gamma}^{(i)} = \Delta\boldsymbol{R}^{(i-1)} - \Delta\boldsymbol{R}^{(i)}
$$
Condition Number $c^{(i)}$ Calculation [*Finite Element Procedures (2nd), P760, (8.104)*].
$$
c^{(i)} = \sqrt{\frac{\big(\boldsymbol{\delta}^{(i)}\big)^{T} \boldsymbol{\gamma}^{(i)}}{\beta \big(\boldsymbol{\delta}^{(i)}\big)^{T}  \Delta\boldsymbol{R}^{(i-1)}}}
$$
Matrix Update Matrix $\boldsymbol{A}^{(i)}$ Calculation [*Finite Element Procedures (2nd), P760, (8.101), (8.102), (8.103)*].
$$
{}^{t + \Delta t}\boldsymbol{K}^{(i-1)} \boldsymbol{\delta}^{(i)} = \beta \Delta\boldsymbol{R}^{(i-1)}
$$

$$
\boldsymbol{v}^{(i)} = 
-\beta \ c^{(i)} \Delta\boldsymbol{R}^{(i-1)} -
\boldsymbol{\gamma}^{(i)}
$$

$$
\boldsymbol{w}^{(i)} = \frac{\boldsymbol{\delta}^{(i)}}{\big(\boldsymbol{\delta}^{(i)}\big)^{T} \boldsymbol{\gamma}^{(i)}}
$$

$$
\boldsymbol{A}^{(i)} = I + \boldsymbol{v}^{(i)} \Big(\boldsymbol{w}^{(i)}\Big)^{T}
$$

Coefficient Matrix $\Big( {}^{t + \Delta t}\boldsymbol{K}^{-1} \Big)^{(i)}$ Update [*Finite Element Procedures (2nd), P760, (8.100)*].
$$
\Big( {}^{t + \Delta t}\boldsymbol{K}^{-1} \Big)^{(i)} =
\begin{cases}
\Big(\boldsymbol{A}^{(i)}\Big)^{T} 
\Big( {}^{t + \Delta t}\boldsymbol{K}^{-1} \Big)^{(i-1)} 
\boldsymbol{A}^{(i)} & c^{(i)} \leq PTOL \\
\Big( {}^{t + \Delta t}\boldsymbol{K}^{-1} \Big)^{(i-1)} & \text{otherwise}
\end{cases}
$$

3. End of Increment

$$
{}^{t + \Delta t}\boldsymbol{F}
$$

$$
{}^{t + \Delta t}\boldsymbol{U}
$$

$$
{}^{\tau}\boldsymbol{K}^{-1}
$$

#### Load-Displacement-Constraint Iterator

1. Beginning of Increment

$$
{}^{t + \Delta t}\boldsymbol{F}^{(0)} = {}^{t}\boldsymbol{F}
$$

$$
{}^{t + \Delta t}\boldsymbol{U}^{(0)} = {}^{t}\boldsymbol{U}
$$

Determine whether to update tangent stiffness matrix [*Finite Element Procedures (2nd), P759*] or use other self-adaptive ways.
$$
{}^{\tau}\boldsymbol{K} =
\frac{\partial \boldsymbol{F}}{\partial \boldsymbol{U}} \Bigg|_{{}^{t} \boldsymbol{U}}
$$
Evaluate ${}^{t + \Delta t}\boldsymbol{R}$.

2. Each Iteration

Solve for $\Delta \boldsymbol{\bar U}^{(i)}$. [*Finite Element Procedures (2nd), P764, (8.113)*]
$$
{}^{\tau}\boldsymbol{K} \Delta\boldsymbol{\bar U}^{(i)} =
{}^{t + \Delta t}\lambda^{(i-1)} \ {}^{t + \Delta t}\boldsymbol{R} - {}^{t + \Delta t}\boldsymbol{F}^{(i-1)}
$$
Solve for $\Delta \boldsymbol{\bar{\bar U}}^{(i)}$. [*Finite Element Procedures (2nd), P764, (8.114)*]
$$
{}^{\tau}\boldsymbol{K} \Delta\boldsymbol{\bar{\bar U}}^{(i)} = {}^{t + \Delta t}\boldsymbol{R}
$$
Solve for $\Delta \lambda^{(i)}$ using constraint equation. 

Calculate ${}^{t + \Delta t}\lambda^{(i)}$ [*Finite Element Procedures (2nd), P764, (8.116)*].
$$
{}^{t + \Delta t}\lambda^{(i)} = {}^{t + \Delta t}\lambda^{(i-1)} + \Delta\lambda^{(i)}
$$
Calculate ${}^{t + \Delta t}\boldsymbol{U}^{(i)}$ [*Finite Element Procedures (2nd), P764, (8.117)*].
$$
{}^{t + \Delta t}\boldsymbol{U}^{(i)} = {}^{t + \Delta t}\boldsymbol{U}^{(i-1)} +
\Delta\boldsymbol{\bar U}^{(i)} +
\Delta\lambda^{(i)} \Delta\boldsymbol{\bar{\bar U}}^{(i)}
$$
Evaluate nodal point force vectors ${}^{t + \Delta t}\boldsymbol{F}^{(i)}$.
$$
{}^{t + \Delta t}\boldsymbol{F}^{(i)}
= \boldsymbol{F} \Big|_{{}^{t + \Delta t} \boldsymbol{U}^{(i)}}
$$

3. End of Increment

Store:
$$
{}^{t + \Delta t}\boldsymbol{F}
$$

$$
{}^{t + \Delta t}\boldsymbol{U}
$$

$$
{}^{\tau}\boldsymbol{K}
$$

##### Constant Increment of External Work Criterion

1. Beginning of Increment

Adaptively select $W$ [*Finite Element Procedures (2nd), P764*].

2. Each Iteration

Solve for $\Delta \lambda^{(i)}$ [*Finite Element Procedures (2nd), P764, (8.112), (8.118)*].
$$
\Delta\lambda^{(i)} = 
\begin{cases}
2W \Big/ \ \Big({}^{t + \Delta t}\boldsymbol{R}^{T} \Delta \boldsymbol{U}^{(i)}\Big) - 2 \ {}^{t}\lambda & i=1 \\
- \boldsymbol{R}^{T} \Delta \boldsymbol{\bar U}^{(i)} \big/ \boldsymbol{R}^{T} \Delta \boldsymbol{\bar{\bar U}}^{(i)} & i=2,3...
\end{cases}
$$

##### Spherical Constant Arc Length Criterion

1. Beginning of Increment

Adaptively select $\Delta l$ [*Finite Element Procedures (2nd), P764*].

2. Each Iteration

[*Finite Element Procedures Solution Manual (2nd), P442, 8.35*].
$$
\boldsymbol{a} = {}^{t + \Delta t}\boldsymbol{U}^{(i-1)} + \Delta\bar{\boldsymbol{U}}^{(i)}
$$
Search for $\beta$ such that $d^{2} > 0$. 
$$
b = \sqrt{\beta + \Big|\Big|\Delta\boldsymbol{\bar{\bar U}}^{(i)}\Big|\Big|_{2}^{2}}
$$

$$
c = \Bigg( \beta \lambda^{(i-1)} + \boldsymbol{a}^{T} \Delta\boldsymbol{\bar{\bar U}}^{(i)} \Bigg) \bigg/ b
$$

$$
d^{2} = \beta \big(\Delta l\big)^{2} - \beta \Big(\lambda^{(i-1)}\Big)^{2} + ||\boldsymbol{a}||_{2}^{2} + c^{2}
$$

Solve for $\Delta \lambda^{(i)}$.
$$
\Delta \lambda^{(i)} = (d-c) / b
$$

### Criterion

#### Displacement Criterion

[*Finite Element Procedures (2nd), P764, (8.119)*].
$$
\Big|\Big| \Delta\boldsymbol{U}^{(i)} \Big|\Big|_{2} \leq 
\epsilon_{D} \Big|\Big| {}^{t + \Delta t}\boldsymbol{U}^{(i)} \Big|\Big|_{2}
$$

#### Unbalanced Force Criterion

[*Finite Element Procedures (2nd), P765, (8.120)*].
$$
\Big|\Big| {}^{t + \Delta t}\boldsymbol{R} - {}^{t + \Delta t}\boldsymbol{F}^{(i)} \Big|\Big|_{2} \leq 
\epsilon_{F} \Big|\Big| {}^{t + \Delta t}\boldsymbol{R} - {}^{t}\boldsymbol{F} \Big|\Big|_{2}
$$

#### Internal Energy Criterion

[*Finite Element Procedures (2nd), P765, (8.121)*].
$$
\Big(\Delta\boldsymbol{U}^{(i)}\Big)^{T} \Big( {}^{t + \Delta t}\boldsymbol{R} - {}^{t + \Delta t}\boldsymbol{F}^{(i)} \Big) \leq
\epsilon_{E} \Big(\Delta\boldsymbol{U}^{(1)}\Big)^{T} \Big( {}^{t + \Delta t}\boldsymbol{R} - {}^{t}\boldsymbol{F} \Big)
$$

### Searcher

Get the improvement Solution ${}^{t + \Delta t}\boldsymbol{U}^{(i)}$. Either not use line search [*Finite Element Procedures (2nd), P755, (8.90)*].
$$
{}^{t + \Delta t}\boldsymbol{U}^{(i)} =
{}^{t + \Delta t} \boldsymbol{U}^{(i-1)} + \Delta \boldsymbol{U}^{(i)}
$$
Or use line search [*Finite Element Procedures (2nd), P760, (8.98)*] to compromise between single iteration expense and total number of iterations for convergence.
$$
{}^{t + \Delta t}\boldsymbol{U}^{(i)} =
{}^{t + \Delta t} \boldsymbol{U}^{(i-1)} + \beta \Delta \boldsymbol{U}^{(i)}
$$
Search for $\beta$. [*Finite Element Procedures (2nd), P760, (8.99)*]
$$
\Big(\Delta\boldsymbol{U}^{(i)}\Big)^{T} \Big( {}^{t + \Delta t}\boldsymbol{R} - {}^{t + \Delta t}\boldsymbol{F}^{(i)} \Big) \leq
STOL \Big(\Delta\boldsymbol{U}^{(i)}\Big)^{T} \Big( {}^{t + \Delta t}\boldsymbol{R} - {}^{t + \Delta t}\boldsymbol{F}^{(i-1)} \Big)
$$

### Integrator

**Interpolation**
$$
u_{i} = \hat h_{n} \hat U_{in}
$$

$$
{}^{t}u_{i} = \hat h_{n} \ {}^{t}\hat U_{in}
$$

$$
x_{i} = \hat h_{n} \hat X_{in}
$$

**Jacobi Matrix**
$$
\boldsymbol{\hat J} = \left(\begin{matrix}
x_{0,0} & x_{0,1} & x_{0,2} \\
x_{1,0} & x_{1,1} & x_{1,2} \\
x_{2,0} & x_{2,1} & x_{2,2}
\end{matrix}\right)
$$

$$
\hat J_{ij} 
= x_{i,j}
= \frac{\partial x_{i}}{\partial r_{j}}
= \hat h_{n,i} \hat X_{in}
$$

$$
\det\boldsymbol{\hat J} 

= \epsilon_{ijk} \ x_{i,0} \ x_{j,1} \ x_{j,2} 

= (x_{1,0} x_{2,1} - x_{1,1} x_{2,0}) x_{0,2} +
(x_{0,1} x_{2,0} - x_{0,0} x_{2,1}) x_{1,2} +
(x_{0,0} x_{1,1} - x_{0,1} x_{1,0}) x_{2,2}
$$

$$
\boldsymbol{\hat J}^{*} = 

\left(\begin{matrix}
x_{1,1} x_{2,2} - x_{1,2} x_{2,1} & x_{0,2} x_{2,1} - x_{0,1} x_{2,2} & x_{0,1} x_{1,2} - x_{0,2} x_{1,1} \\
x_{1,2} x_{2,0} - x_{1,0} x_{2,2} & x_{0,0} x_{2,2} - x_{0,2} x_{2,0} & x_{0,2} x_{1,0} - x_{0,0} x_{1,2} \\
x_{1,0} x_{2,1} - x_{1,1} x_{2,0} & x_{0,1} x_{2,0} - x_{0,0} x_{2,1} & x_{0,0} x_{1,1} - x_{0,1} x_{1,0}
\end{matrix}\right)
$$

$$
\boldsymbol{\hat J}^{-1} = 
\frac{1}{\det \boldsymbol{\hat J}} \boldsymbol{\hat J^{*}}
$$

**Derivative**
$$
{}_{0}\hat h_{n,j}
= \frac{\partial \hat h_{n}}{\partial \ {}^{0}x_{j}}
= {}_{0}\hat J_{jl}^{-1} \frac{\partial \hat h_{n}}{\partial r_{l}}
= {}_{0}\hat J_{jl}^{-1} \hat h_{n,l}
$$

$$
{}_{0}u_{i,j}
= \frac{\partial u_{i}}{\partial \ {}^{0}x_{j}}
= {}_{0}\hat h_{n,j} \hat U_{in}
= {}_{0}\hat J_{jl}^{-1} \hat h_{n,l} \hat U_{in}
$$

$$
{}_{0}^{t}u_{i,j}
= \frac{\partial \ {}^{t}u_{i}}{\partial \ {}^{0}x_{j}}
= {}_{0}\hat h_{n,j} \ {}^{t}\hat U_{in}
= {}_{0}\hat J_{jl}^{-1} \hat h_{n,l} \ {}^{t}\hat U_{in}
$$

**Strain Increment**
$$
{}_{0}\eta_{ij} 
= \frac{1}{2} {}_{0}u_{k,i} \ {}_{0}u_{k,j}
= \frac{1}{2} \ {}^{0}\hat J_{il}^{-1} \hat h_{n,l} \ {}^{0}\hat J_{jl'}^{-1} \hat h_{n',l'} \hat U_{kn} \hat U_{kn'}
$$

$$
\begin{align}

{}_{0}e_{ij}

&= \frac{1}{2} \Big( {}_{0}u_{i,j} + {}_{0}u_{j,i} + {}_{0}^{t}u_{k,i} \ {}_{0}u_{k,j} + {}_{0}u_{k,i} \ {}_{0}^{t}u_{k,j} \Big) \\

&= \frac{1}{2} \Big( {}^{0}\hat J_{jl}^{-1} \hat h_{n,l} \hat U_{in} +
{}^{0}\hat J_{il}^{-1} \hat h_{n,l} \hat U_{jn} + 
{}^{0}\hat J_{il}^{-1} \hat h_{n,l} \ {}^{0}\hat J_{jl'}^{-1} \hat h_{n',l'} \ {}^{t}\hat U_{kn} \ \hat U_{kn'} +
{}^{0}\hat J_{il}^{-1} \hat h_{n,l} \ {}^{0}\hat J_{jl'}^{-1} \hat h_{n',l'} \ \hat U_{kn} \ {}^{t}\hat U_{kn'} \Big)

\end{align}
$$

**Virtual Strain Increment**
$$
\begin{align}

\delta \ {}_{0}\eta_{ij}

&= {}_{0}\eta_{ij,pm} \delta \hat U_{pm} \\

&= \frac{1}{2} \ {}^{0}\hat J_{il}^{-1} \hat h_{n,l} \ {}^{0}\hat J_{jl'}^{-1} \hat h_{n',l'} \frac{\partial \big(\hat U_{kn} \hat U_{kn'}\big)}{\partial \hat U_{pm}} \delta \hat U_{pm} \\

&= \frac{1}{2} \ {}^{0}\hat J_{il}^{-1} \hat h_{n,l} \ {}^{0}\hat J_{jl'}^{-1} \hat h_{n',l'}
\Bigg(
\delta_{n'm} \delta_{pk} \hat U_{kn}
+ \delta_{nm} \delta_{pk} \hat U_{kn'}
\Bigg) \delta \hat U_{pm} \\

&= \frac{1}{2} \ {}^{0}\hat J_{il}^{-1} \hat h_{n,l} \ {}^{0}\hat J_{jl'}^{-1} \hat h_{n',l'}
\Bigg(
\delta_{n'm} \hat U_{pn}
+ \delta_{nm} \hat U_{pn'}
\Bigg) \delta \hat U_{pm} \\

&= \frac{1}{2} \Bigg (
{}^{0}\hat J_{il}^{-1} \hat h_{n,l} \ {}^{0}\hat J_{jl'}^{-1} \hat h_{m,l'} \hat U_{pn}
+ {}^{0}\hat J_{il}^{-1} \hat h_{m,l} \ {}^{0}\hat J_{jl'}^{-1} \hat h_{n',l'} \hat U_{pn'}
\Bigg) \delta \hat U_{pm} \\

&= \frac{1}{2} \Bigg (
{}^{0}\hat J_{il}^{-1} \hat h_{n,l} \ {}^{0}\hat J_{jl'}^{-1} \hat h_{m,l'} \hat U_{pn}
+ {}^{0}\hat J_{il}^{-1} \hat h_{m,l} \ {}^{0}\hat J_{jl'}^{-1} \hat h_{n,l'} \hat U_{pn}
\Bigg) \delta \hat U_{pm} \\

&= \frac{1}{2} \Bigg( {}^{0}\hat J_{il}^{-1} \hat h_{n,l} \ {}^{0}\hat J_{jl'}^{-1} \hat h_{m,l'} +
{}^{0}\hat J_{il}^{-1} \hat h_{m,l} \ {}^{0}\hat J_{jl'}^{-1} \hat h_{n,l'} \Bigg)
\hat U_{pn} \delta \hat U_{pm}

\end{align}
$$

$$
\begin{align}

\delta \ {}_{0}e_{ij}

&= {}_{0}e_{ij,pm} \delta \hat U_{pm} \\

&= \frac{1}{2} \Big( {}^{0}\hat J_{jl}^{-1} \hat h_{n,l} \frac{\partial \hat U_{in}}{\partial \hat U_{pm}} +
{}^{0}\hat J_{il}^{-1} \hat h_{n,l} \frac{\partial \hat U_{jn}}{\partial \hat U_{pm}} + 
{}^{0}\hat J_{il}^{-1} \hat h_{n,l} \ {}^{0}\hat J_{jl'}^{-1} \hat h_{n',l'} \frac{\partial \hat U_{kn'}}{\partial \hat U_{pm}} {}^{t}\hat U_{kn} +
{}^{0}\hat J_{il}^{-1} \hat h_{n,l} \ {}^{0}\hat J_{jl'}^{-1} \hat h_{n',l'} \frac{\partial \hat U_{kn}}{\partial \hat U_{pm}} {}^{t}\hat U_{kn'} \Bigg) \delta \hat U_{pm} \\

&= \frac{1}{2} \Big( {}^{0}\hat J_{jl}^{-1} \hat h_{n,l} \delta_{ip} \delta_{nm} +
{}^{0}\hat J_{il}^{-1} \hat h_{n,l} \delta_{jp} \delta_{nm} + 
{}^{0}\hat J_{il}^{-1} \hat h_{n,l} \ {}^{0}\hat J_{jl'}^{-1} \hat h_{n',l'} \delta_{kp} \delta_{n'm} \ {}^{t}\hat U_{kn} +
{}^{0}\hat J_{il}^{-1} \hat h_{n,l} \ {}^{0}\hat J_{jl'}^{-1} \hat h_{n',l'} \delta_{kp} \delta_{nm} \ {}^{t}\hat U_{kn'} \Bigg) \delta \hat U_{pm} \\

&= \frac{1}{2} \Big( {}^{0}\hat J_{jl}^{-1} \hat h_{m,l} \delta_{ip} +
{}^{0}\hat J_{il}^{-1} \hat h_{m,l} \delta_{jp} + 
{}^{0}\hat J_{il}^{-1} \hat h_{n,l} \ {}^{0}\hat J_{jl'}^{-1} \hat h_{m,l'} \ {}^{t}\hat U_{pn} +
{}^{0}\hat J_{il}^{-1} \hat h_{m,l} \ {}^{0}\hat J_{jl'}^{-1} \hat h_{n',l'} \ {}^{t}\hat U_{pn'} \Bigg) \delta \hat U_{pm} \\

&= \frac{1}{2} \Big( {}^{0}\hat J_{jl}^{-1} \hat h_{m,l} \delta_{ip} +
{}^{0}\hat J_{il}^{-1} \hat h_{m,l} \delta_{jp} + 
{}^{0}\hat J_{il}^{-1} \hat h_{n,l} \ {}^{0}\hat J_{jl'}^{-1} \hat h_{m,l'} {}^{t}\hat U_{pn} +
{}^{0}\hat J_{il}^{-1} \hat h_{m,l} \ {}^{0}\hat J_{jl'}^{-1} \hat h_{n,l'} {}^{t}\hat U_{pn} \Bigg) \delta \hat U_{pm}

\end{align}
$$

**Virtual Displacement Increment**
$$
\delta u_{i}

= u_{i,pm} \ \delta \hat U_{pm}

= \hat h_{n} \frac{\partial \hat U_{in}}{\partial \hat U_{pm}} \delta \hat U_{pm}

= \hat h_{n} \delta_{ip} \delta_{nm} \delta \hat U_{pm}

= \hat h_{m} \delta_{} \hat U_{im}
$$
**Volume Element**

[*https://en.wikipedia.org/wiki/Volume_element*]
$$
\text{d} \ {}^{0}\hat V = \Big|\det \boldsymbol{{}^{0}\hat J}\Big| \ \text{d}r_{1} \ \text{d}r_{2} \ \text{d}r_{3}
$$
**Global Degree Assembly**
$$
\hat U_{pm} = U_{\hat f(m,p)}
$$

$$
\delta \hat U_{pm}
= \frac{\partial \hat U_{pm}}{\partial U_{q}} \delta U_{q}
= \delta_{\hat f(m,p)q} \delta U_{q}
= \delta U_{\hat f(m,p)}
$$

**U**

#### Nodal Point Force Vectors Integrator

$$
\iiint_{{}^{0}V}
{}_{0}^{t}S_{ij} \ \delta \ {}_{0}e_{ij} \ \text{d} \ {}^{0}V =

\sum_{e=1}^{E} \iiint_{{}^{0}\hat V}
{}_{0}^{t}S_{ij} \delta \ {}_{0}e_{ij} \text{d} \ {}^{0}\hat V
$$

Since ${}_{0}^{t} S_{ij} = {}_{0}^{t} S_{ji}$
$$
\begin{align}

{}_{0}^{t}S_{ij} \delta \ {}_{0}e_{ij}

&= \frac{1}{2} \Big( {}_{0}^{t}S_{ij} \ {}^{0}\hat J_{jl}^{-1} \hat h_{m,l} \delta_{ip} +
{}_{0}^{t}S_{ij} \ {}^{0}\hat J_{il}^{-1} \hat h_{m,l} \delta_{jp} +
{}_{0}^{t}S_{ij} \ {}^{0}\hat J_{il}^{-1} \hat h_{n,l} \ {}^{0}\hat J_{jl'}^{-1} \hat h_{m,l'} {}^{t}\hat U_{pn} +
{}_{0}^{t}S_{ij} \ {}^{0}\hat J_{il}^{-1} \hat h_{m,l} \ {}^{0}\hat J_{jl'}^{-1} \hat h_{n,l'} {}^{t}\hat U_{pn} \Bigg) \delta \hat U_{pm} \\

&= {}_{0}^{t}S_{ij} 
\Big( {}^{0}\hat J_{jl}^{-1} \hat h_{m,l} \delta_{ip} +
{}^{0}\hat J_{il}^{-1} \hat h_{n,l} \ {}^{0}\hat J_{jl'}^{-1} \hat h_{m,l'} {}^{t}\hat U_{pn} \Big)
\delta \hat U_{pm}

\end{align}
$$
Thus
$$
\iiint_{{}^{0}\hat V} {}_{0}^{t}S_{ij} \delta \ {}_{0}e_{ij} \text{d} \ {}^{0}\hat V

= \delta \hat U_{pm} \iiint
{}_{0}^{t}S_{ij} 
\Big( {}^{0}\hat J_{jl}^{-1} \hat h_{m,l} \delta_{ip} +
{}^{0}\hat J_{il}^{-1} \hat h_{n,l} \ {}^{0}\hat J_{jl'}^{-1} \hat h_{m,l'} {}^{t}\hat U_{pn} \Big) \Big|\det \boldsymbol{{}^{0}\hat J}\Big| \ \text{d}r_{1} \ \text{d}r_{2} \ \text{d}r_{3} \
$$


#### Tangent Stiffness Matrices Integrator

##### Linear Strain Incremental Stiffness Matrices Integrator

[*Finite Element Procedures (2nd), P524, TABLE 6.2*]
$$
\iiint_{{}^{0}V} {}_{0}C_{ijrs} \ {}_{0}e_{rs} \ \delta \ {}_{0}e_{ij} \ \text{d} \ {}^{0}V =
\sum_{e=1}^{E} \iiint_{{}^{0}\hat V} {}_{0}C_{ijrs} \ {}_{0}e_{rs} \ \delta \ {}_{0}e_{ij} \ \text{d} \ {}^{0}\hat V
$$

Since ${}_{0}C_{ijrs} = {}_{0}C_{ijsr}$
$$
\begin{align}

{}_{0}C_{ijrs} \ {}_{0}e_{rs}

=& \frac{1}{2} \Big( {}_{0}C_{ijrs} \ {}^{0}\hat J_{sl}^{-1} \hat h_{n,l} \hat U_{rn} +
{}_{0}C_{ijrs} \ {}^{0}\hat J_{rl}^{-1} \hat h_{n,l} \hat U_{sn} + \\
& {}_{0}C_{ijrs} \ {}^{0}\hat J_{rl}^{-1} \hat h_{n,l} \ {}^{0}\hat J_{sl'}^{-1} \hat h_{n',l'} \ {}^{t}\hat U_{kn} \ \hat U_{kn'} +
{}_{0}C_{ijrs} \ {}^{0}\hat J_{rl}^{-1} \hat h_{n,l} \ {}^{0}\hat J_{sl'}^{-1} \hat h_{n',l'} \ \hat U_{kn} \ {}^{t}\hat U_{kn'} \Big) \\

=& {}_{0}C_{ijrs} \ {}^{0}\hat J_{sl}^{-1} \hat h_{n,l} \hat U_{rn} +
{}_{0}C_{ijrs} \ {}^{0}\hat J_{rl}^{-1} \hat h_{n,l} \ {}^{0}\hat J_{sl'}^{-1} \hat h_{n',l'} \ {}^{t}\hat U_{kn} \ \hat U_{kn'} \\

=& {}_{0}C_{ijrs} \Big( {}^{0}\hat J_{sl}^{-1} \hat h_{n,l} \hat U_{rn} +
{}^{0}\hat J_{rl}^{-1} \hat h_{n,l} \ {}^{0}\hat J_{sl'}^{-1} \hat h_{n',l'} \ {}^{t}\hat U_{kn} \ \hat U_{kn'} \Big)

\end{align}
$$
Since ${}_{0}C_{ijrs} = {}_{0}C_{jirs}$
$$
\begin{align}

{}_{0}C_{ijrs} \ {}_{0}e_{rs} \ \delta \ {}_{0}e_{ij}

=&  \frac{1}{2} {}_{0}C_{ijrs} \Big( {}^{0}\hat J_{sl}^{-1} \hat h_{n,l} \hat U_{rn} +
{}^{0}\hat J_{rl}^{-1} \hat h_{n,l} \ {}^{0}\hat J_{sl'}^{-1} \hat h_{n',l'} \ {}^{t}\hat U_{kn} \ \hat U_{kn'} \Big) \\
& \Big( {}^{0}\hat J_{jl''}^{-1} \hat h_{m,l''} \delta_{ip} +
{}^{0}\hat J_{il''}^{-1} \hat h_{m,l''} \delta_{jp} + 
{}^{0}\hat J_{il''}^{-1} \hat h_{n'',l''} \ {}^{0}\hat J_{jl'''}^{-1} \hat h_{m,l'''} {}^{t}\hat U_{pn''} +
{}^{0}\hat J_{il''}^{-1} \hat h_{m,l''} \ {}^{0}\hat J_{jl'''}^{-1} \hat h_{n'',l'''} {}^{t}\hat U_{pn''} \Bigg) \delta \hat U_{pm} \\

=& \frac{1}{2} \Big( {}^{0}\hat J_{sl}^{-1} \hat h_{n,l} \hat U_{rn} +
{}^{0}\hat J_{rl}^{-1} \hat h_{n,l} \ {}^{0}\hat J_{sl'}^{-1} \hat h_{n',l'} \ {}^{t}\hat U_{kn} \ \hat U_{kn'} \Big) \\
& \Big( {}_{0}C_{ijrs} \ {}^{0}\hat J_{jl''}^{-1} \hat h_{m,l''} \delta_{ip} +
{}_{0}C_{ijrs} \ {}^{0}\hat J_{il''}^{-1} \hat h_{m,l''} \delta_{jp} + \\
& {}_{0}C_{ijrs} \ {}^{0}\hat J_{il''}^{-1} \hat h_{n'',l''} \ {}^{0}\hat J_{jl'''}^{-1} \hat h_{m,l'''} {}^{t}\hat U_{pn''} +
{}_{0}C_{ijrs} \ {}^{0}\hat J_{il''}^{-1} \hat h_{m,l''} \ {}^{0}\hat J_{jl'''}^{-1} \hat h_{n'',l'''} {}^{t}\hat U_{pn''} \Bigg) \delta \hat U_{pm} \\

=& \Big( {}^{0}\hat J_{sl}^{-1} \hat h_{n,l} \hat U_{rn} +
{}^{0}\hat J_{rl}^{-1} \hat h_{n',l} \ {}^{0}\hat J_{sl'}^{-1} \hat h_{n,l'} \ {}^{t}\hat U_{kn'} \ \hat U_{kn} \Big) \\
& \Big( {}_{0}C_{ijrs} \ {}^{0}\hat J_{jl''}^{-1} \hat h_{m,l''} \delta_{ip} + 
{}_{0}C_{ijrs} \ {}^{0}\hat J_{il''}^{-1} \hat h_{n'',l''} \ {}^{0}\hat J_{jl'''}^{-1} \hat h_{m,l'''} {}^{t}\hat U_{pn''} \Big) \delta \hat U_{pm} \\

=& \Big( {}^{0}\hat J_{sl}^{-1} \hat h_{n,l} \delta_{kr} +
{}^{0}\hat J_{rl}^{-1} \hat h_{n',l} \ {}^{0}\hat J_{sl'}^{-1} \hat h_{n,l'} \ {}^{t}\hat U_{kn'} \Big)
{}_{0}C_{ijrs} \Big( {}^{0}\hat J_{jl''}^{-1} \hat h_{m,l''} \delta_{ip} + 
{}^{0}\hat J_{il''}^{-1} \hat h_{n'',l''} \ {}^{0}\hat J_{jl'''}^{-1} \hat h_{m,l'''} {}^{t}\hat U_{pn''} \Big)
\hat U_{kn} \ \delta \hat U_{pm} 

\end{align}
$$
Thus
$$
\begin{align}
\iiint_{{}^{0}\hat V} {}_{0}C_{ijrs} \ {}_{0}e_{rs} \ \delta \ {}_{0}e_{ij}

=& \hat U_{kn} \ \delta \hat U_{pm} \iiint
\Big( {}^{0}\hat J_{sl}^{-1} \hat h_{n,l} \delta_{kr} +
{}^{0}\hat J_{rl}^{-1} \hat h_{n',l} \ {}^{0}\hat J_{sl'}^{-1} \hat h_{n,l'} \ {}^{t}\hat U_{kn'} \Big)
{}_{0}C_{ijrs} \\
&\Big( {}^{0}\hat J_{jl''}^{-1} \hat h_{m,l''} \delta_{ip} + 
{}^{0}\hat J_{il''}^{-1} \hat h_{n'',l''} \ {}^{0}\hat J_{jl'''}^{-1} \hat h_{m,l'''} {}^{t}\hat U_{pn''} \Big)
\Big|\det \boldsymbol{{}^{0}\hat J}\Big| \ \text{d}r_{1} \ \text{d}r_{2} \ \text{d}r_{3}

\end{align}
$$

##### Nonlinear Strain Incremental Stiffness Matrices Integrator

[*Finite Element Procedures (2nd), P524, TABLE 6.2*]
$$
\iiint_{{}^{0}V} {}_{0}^{t}S_{ij} \ \delta \ {}_{0}\eta_{ij} \ \text{d} \ {}^{0}V =

\sum_{e=1}^{E} \iiint_{{}^{0}\hat V} {}_{0}^{t}S_{ij} \ \delta \ {}_{0}\eta_{ij} \ \text{d} \ {}^{0}\hat V
$$

Since ${}_{0}^{t} S_{ij} = {}_{0}^{t} S_{ji}$
$$
\begin{align}

{}_{0}^{t}S_{ij} \ \delta \ {}_{0} \eta_{ij}

&= \frac{1}{2} \Bigg(
{}_{0}^{t}S_{ij} \ {}^{0}\hat J_{il}^{-1} \hat h_{n,l} \ {}^{0}\hat J_{jl'}^{-1} \hat h_{m,l'}
+ {}_{0}^{t}S_{ij} \ {}^{0}\hat J_{il}^{-1} \hat h_{m,l} \ {}^{0}\hat J_{jl'}^{-1} \hat h_{n,l'}
\Bigg) \hat U_{pn} \ \delta \hat U_{pm} \\

&= \frac{1}{2} \Bigg(
{}_{0}^{t}S_{ij} \ {}^{0}\hat J_{il}^{-1} \hat h_{n,l} \ {}^{0}\hat J_{jl'}^{-1} \hat h_{m,l'}
+ {}_{0}^{t}S_{ji} \ {}^{0}\hat J_{jl'}^{-1} \hat h_{m,l'} \ {}^{0}\hat J_{il}^{-1} \hat h_{n,l}
\Bigg) \hat U_{pn} \ \delta \hat U_{pm} \\

&= {}_{0}^{t}S_{ij} \ {}^{0}\hat J_{il}^{-1} \hat h_{n,l} \ {}^{0}\hat J_{jl'}^{-1} \hat h_{m,l'} \ \hat U_{pn} \ \delta \hat U_{pm}

\end{align}
$$

Thus

$$
\iiint_{{}^{0} \hat V} {}_{0}^{t}S_{ij} \ \delta \ {}_{0} \eta_{ij} \ \text{d} \ {}^{0}\hat V 

= \hat U_{pn} \ \delta \hat U_{pm}
\iiint
{}_{0}^{t}S_{ij} \ {}^{0}\hat J_{il}^{-1} \hat h_{n,l} \ {}^{0}\hat J_{jl'}^{-1} \hat h_{m,l'} \
\Big|\det\boldsymbol{{}^{0}\hat J}\Big| \ \text{d}r_{1} \ \text{d}r_{2} \ \text{d}r_{3} \
$$

#### External Nodal Load Vectors Integrator

##### Displacement-independent Body Load Vectors

[*Finite Element Procedures (2nd), P542, TABLE 6.4*]
$$
\iiint_{{}^{t+\Delta t}V} {}_{0}^{t + \Delta t}f_{i}^{B} \ \delta u_{i} \ \text{d} {}^{0}\hat V

= \sum_{e=1}^{E} \iiint_{{}^{0}\hat V} {}_{0}^{t + \Delta t}f_{i}^{B} \ \delta u_{i} \ \text{d} {}^{0}\hat V
$$

$$
{}_{0}^{t + \Delta t}f_{i}^{B} \ \delta u_{i}

= {}_{0}^{t + \Delta t}f_{i}^{B} \ \hat h_{m} \delta \hat U_{im}
$$

Thus
$$
\begin{align}

\iiint_{{}^{0}\hat V} {}_{0}^{t + \Delta t}f_{i}^{B} \ \delta u_{i} \ \text{d} {}^{0}\hat V

= \delta \hat U_{im} \iiint
{}_{0}^{t + \Delta t}f_{i}^{B} \ \hat h_{m} \
\Big|\det\boldsymbol{{}^{0}\hat J}\Big| \ \text{d}r_{1} \ \text{d}r_{2} \ \text{d}r_{3}
 
\end{align}
$$

##### Space Attached Displacement-dependent Surface Load Vectors

[*Finite Element Procedures (2nd), P528, (6.84)*]
$$
\iint_{{}^{t + \Delta t}S_{f}} {}^{t + \Delta t}f_{i}^{S} \ \delta u_{i}^{S} \ \text{d} \ {}^{t + \Delta t}S

= \sum_{e=1}^{E} \iint_{{}^{t + \Delta t}\hat S_{f}} {}^{t + \Delta t}f_{i}^{S} \ \delta u_{i}^{S} \ \text{d} \ {}^{t + \Delta t}\hat S
$$

Interpolation.
$$
{}^{t}x_{i}^{S} = \hat g_{n} \ {}^{t}\hat X_{in}^{S}
$$

$$
u_{i}^{S} = \hat g_{n} \ \hat U_{in}^{S}
$$

Derivative.
$$
{}^{t}x_{i,q}^{S}

= \frac{\partial \ {}^{t}x_{i}^{S}}{\partial s_{q}}

= \frac{\partial \hat g_{n}}{\partial s_{q}} \ {}^{t}\hat X_{in}^{S}

= \hat g_{n,q} \ {}^{t}\hat X_{in}^{S}
$$

$$
u_{i,q}^{S}

= \frac{\partial u_{i}^{S}}{\partial s_{q}}

= \frac{\partial \hat g_{n}}{\partial s_{q}} \hat U_{in}^{S}

= \hat g_{n,q} \hat U_{in}^{S}
$$

Virtual Displacement Increment.
$$
\delta u_{i}^{S}

= u_{i,pm}^{S} \ \delta \hat U_{pm}^{S}

= \hat g_{n} \frac{\partial \hat U_{in}^{S}}{\partial \hat U_{pm}^{S}} \delta \hat U_{pm}^{S}

= \hat g_{n} \delta_{ip} \delta_{nm} \delta \hat U_{pm}^{S}

= \hat g_{m} \delta_{} \hat U_{im}^{S}
$$
Product of Normal and Surface Element, omit nonlinear terms [*Schwlizerhof, Karl and E. Ramm. “Displacement dependent pressure loads in nonlinear finite element analyses.” (1984). (18)*].

[*https://en.wikipedia.org/wiki/Surface_integral*]
$$
\begin{align}
{}^{t + \Delta t}n_{i} \ \text{d} \ {}^{t+\Delta t}\hat S

&= e_{ijk} \frac{\partial \ {}^{t + \Delta t}x_{j}^{S}}{\partial s_{1}} \frac{\partial \ {}^{t + \Delta t}x_{k}^{S}}{\partial s_{2}} \text{d}s_{1} \ \text{d}s_{2} \\

&= e_{ijk} \frac{\partial \big({}^{t}x_{j}^{S} + u_{j}^{S}\big)}{\partial s_{1}} \frac{\partial \big({}^{t}x_{k}^{S} + u_{k}^{S}\big)}{\partial s_{2}} \text{d}s_{1} \ \text{d}s_{2} \\

&\dot{=} \ e_{ijk} \frac{\partial \ {}^{t}x_{j}^{S}}{\partial s_{1}} \frac{\partial \ {}^{t}x_{k}^{S}}{\partial s_{2}} \text{d}s_{1} \ \text{d}s_{2} +
e_{ijk} \frac{\partial u_{j}^{S}}{\partial s_{1}} \frac{\partial \ {}^{t}x_{k}^{S}}{\partial s_{2}} \text{d}s_{1} \ \text{d}s_{2} +
e_{ijk} \frac{\partial \ {}^{t}x_{j}^{S}}{\partial s_{1}} \frac{\partial u_{k}^{S}}{\partial s_{2}} \text{d}s_{1} \ \text{d}s_{2} \\

&= e_{ijk} \ {}^{t}\hat X_{jn}^{S} \ {}^{t}\hat X_{kn'}^{S} \ \hat g_{n,1} \ \hat g_{n',2} \ \text{d}s_{1} \ \text{d}s_{2} +
e_{ijk} \ \hat U_{jn}^{S} \ {}^{t}\hat X_{kn'}^{S} \ \hat g_{n,1} \ \hat g_{n',2} \ \text{d}s_{1} \ \text{d}s_{2} +
e_{ijk} \ {}^{t}\hat X_{jn}^{S} \ \hat U_{kn'}^{S} \ \hat g_{n,1} \ \hat g_{n',2} \ \text{d}s_{1} \ \text{d}s_{2} \\

&= e_{ijk} \ {}^{t}\hat X_{jn}^{S} \ {}^{t}\hat X_{kn'}^{S} \ \hat g_{n,1} \ \hat g_{n',2} \ \text{d}s_{1} \ \text{d}s_{2} +
e_{ijk} \ \hat U_{jn}^{S} \ {}^{t}\hat X_{kn'}^{S} \Big( \hat g_{n,1} \ \hat g_{n',2} - \hat g_{n',1} \ \hat g_{n,2} \Big) \ \text{d}s_{1} \ \text{d}s_{2}

\end{align}
$$
Surface Traction Component [*Schwlizerhof, Karl and E. Ramm. “Displacement dependent pressure loads in nonlinear finite element analyses.” (1984). (14), Fig. 3.*].
$$
{}^{t + \Delta t}f_{i}^{S}

= {}^{t + \Delta t}\lambda \ {}^{t + \Delta t}f \ {}^{t + \Delta t}n_{i}
$$
Space Attached Surface Traction Magnitude [*Schwlizerhof, Karl and E. Ramm. “Displacement dependent pressure loads in nonlinear finite element analyses.” (1984). (21)*].
$$
{}^{t + \Delta t}f \ 

\dot{=} \ {}^{t}f + \frac{\partial \ {}^{t}f}{\partial \ {}^{t}x_{i}^{S}}u_{i}^{S}

= {}^{t}f + {}_{t}^{t}f_{,i} \ \hat g_{n} \ \hat U_{in}^{S}
$$
External virtual work. [*Schwlizerhof, Karl and E. Ramm. “Displacement dependent pressure loads in nonlinear finite element analyses.” (1984). (17)*].
$$
\begin{align}

\iint_{{}^{t + \Delta t}S_{f}}
{}^{t + \Delta t}f_{i}^{S} \ \delta u_{i}^{S} \ \text{d}\ {}^{t + \Delta t}\hat S

=& \iint_{{}^{t + \Delta t}S_{f}}
{}^{t + \Delta t}\lambda \ {}^{l}f \ \delta u_{i}^{S} \ {}^{t + \Delta t}n_{i} \ \text{d} \ {}^{t + \Delta t}\hat S \\

\\

=& \iint {}^{t + \Delta t}\lambda \ \Bigg( {}^{t}f + {}_{t}^{t}f_{,l} \ \hat g_{n} \ \hat U_{ln}^{S} \Bigg)
\hat g_{m} \delta \hat U_{im}^{S} \\
& \Bigg[
e_{ijk} \ {}^{t}\hat X_{jn}^{S} \ {}^{t}\hat X_{kn'}^{S} \ \hat g_{n,1} \ \hat g_{n',2} \ \text{d}s_{1} \ \text{d}s_{2} +
e_{ijk} \ \hat U_{jn}^{S} \ {}^{t}\hat X_{kn'}^{S} \Big( \hat g_{n,1} \ \hat g_{n',2} - \hat g_{n',1} \ \hat g_{n,2} \Big) \ \text{d}s_{1} \ \text{d}s_{2}
\Bigg] \\

\\

\dot{=} \ & {}^{t + \Delta t}\lambda \ e_{ijk} \ {}^{t}\hat X_{jn}^{S} \ {}^{t}\hat X_{kn'}^{S}\ \delta \hat U_{im}^{S} \ \iint {}^{t}f \ \hat g_{m} \ \hat g_{n,1} \ \hat g_{n',2} \ \text{d}s_{1} \ \text{d}s_{2} + \\

& {}^{t + \Delta t}\lambda \ e_{ijk} \ {}^{t}\hat X_{kn'}^{S}\ \hat U_{jn}^{S} \ \delta \hat U_{im}^{S} \ \iint {}^{t}f \ \hat g_{m} \Big( \hat g_{n,1} \ \hat g_{n',2} - \hat g_{n',1} \ \hat g_{n,2} \Big) \text{d}s_{1} \ \text{d}s_{2} + \\

& {}^{t + \Delta t}\lambda \ e_{ijk} \ {}^{t}\hat X_{jn}^{S} \ {}^{t}\hat X_{kn'}^{S} \ \hat U_{ln}^{S} \ \delta \hat U_{im}^{S} \ \iint {}_{t}^{t}f_{,l} \ \hat g_{m} \ \hat g_{n} \ \hat g_{n,1} \ \hat g_{n',2} \ \text{d}s_{1} \ \text{d}s_{2}

\end{align}
$$

Since
$$
\begin{align}

\iint {}^{t}f \ \hat g_{m} \Big( \hat g_{n,1} \ \hat g_{n',2} - \hat g_{n',1} \ \hat g_{n,2} \Big) \text{d}s_{1} \ \text{d}s_{2}

=& \iint {}^{t}f \ \hat g_{m} \ \hat g_{n,1} \ \hat g_{n',2} \ \text{d}s_{1} \ \text{d}s_{2}
- \iint {}^{t}f \ \hat g_{m} \ \hat g_{n',1} \ \hat g_{n,2} \ \text{d}s_{1} \ \text{d}s_{2}
\\

\\

=& \Bigg(
\int {}^{t}f \ \hat g_{m} \ \hat g_{n} \ \hat g_{n',2} \ \text{d}s_{2} -
\iint {}_{t}^{t}f_{,l} \ {}^{t}x_{l,1}^{S} \ \hat g_{m} \ \hat g_{n} \ \hat g_{n',2} \ \text{d}s_{1} \ \text{d}s_{2} \\
& - \iint {}^{t}f \ \hat g_{m,1} \ \hat g_{n} \ \hat g_{n',2} \ \text{d}s_{1} \ \text{d}s_{2}
- \iint {}^{t}f \ \hat g_{m} \ \hat g_{n} \ \hat g_{n',12} \ \text{d}s_{1} \ \text{d}s_{2}
\Bigg) \\
& - \Bigg(
\int {}^{t}f \ \hat g_{m} \ \hat g_{n} \ \hat g_{n',1} \ \text{d}s_{1} -
\iint {}_{t}^{t}f_{,l} \ {}^{t}x_{l,2}^{S} \ \hat g_{m} \ \hat g_{n',1} \ \hat g_{n} \ \text{d}s_{1} \ \text{d}s_{2} \\
& - \iint {}^{t}f \ \hat g_{m,2} \ \hat g_{n',1} \ \hat g_{n} \ \text{d}s_{1} \ \text{d}s_{2}
- \iint {}^{t}f \ \hat g_{m} \ \hat g_{n',12} \ \hat g_{n} \ \text{d}s_{1} \ \text{d}s_{2}
\Bigg) \\

\\

=& \int {}^{t}f \ \hat g_{m} \ \hat g_{n} \ \hat g_{n',2} \ \text{d}s_{2}
- \int {}^{t}f \ \hat g_{m} \ \hat g_{n} \ \hat g_{n',1} \ \text{d}s_{1} \\
& + \iint {}^{t}f \ \hat g_{n} \ \Big( \hat g_{m,2} \ \hat g_{n',1} - \hat g_{m,1} \ \hat g_{n',2} \Big) \ \text{d}s_{1} \ \text{d}s_{2} \\
& + \iint {}_{t}^{t}f_{,l} \ \hat g_{m} \ \hat g_{n} \ \Big( {}^{t}x_{l,2}^{S} \ \hat g_{n',1} - {}^{t}x_{l,1}^{S} \ \hat g_{n',2} \Big) \ \text{d}s_{1} \ \text{d}s_{2}

\end{align}
$$

Thus
$$
\begin{align}

\iint {}^{t}f \ \hat g_{m} \Big( \hat g_{n,1} \ \hat g_{n',2} - \hat g_{n',1} \ \hat g_{n,2} \Big) \text{d}s_{1} \ \text{d}s_{2}

=& \frac{1}{2} \Bigg[
\iint {}^{t}f \ \hat g_{m} \Big( \hat g_{n,1} \ \hat g_{n',2} - \hat g_{n',1} \ \hat g_{n,2} \Big) \text{d}s_{1} \ \text{d}s_{2} 
+ \iint {}^{t}f \ \hat g_{n} \ \Big( \hat g_{m,2} \ \hat g_{n',1} - \hat g_{m,1} \ \hat g_{n',2} \Big) \ \text{d}s_{1} \ \text{d}s_{2} \\
& + \iint {}_{t}^{t}f_{,l} \ \hat g_{m} \ \hat g_{n} \ \Big( {}^{t}x_{l,2}^{S} \ \hat g_{n',1} - {}^{t}x_{l,1}^{S} \ \hat g_{n',2} \Big) \ \text{d}s_{1} \ \text{d}s_{2} \\
& + \int {}^{t}f \ \hat g_{m} \ \hat g_{n} \ \hat g_{n',2} \ \text{d}s_{2}
- \int {}^{t}f \ \hat g_{m} \ \hat g_{n} \ \hat g_{n',1} \ \text{d}s_{1}
\Bigg] \\

\\

=& \frac{1}{2} \Bigg[
\iint {}^{t}f \ \hat g_{m} \Big( \hat g_{n,1} \ \hat g_{n',2} - \hat g_{n',1} \ \hat g_{n,2} \Big) \text{d}s_{1} \ \text{d}s_{2} 
+ \iint {}^{t}f \ \hat g_{n} \ \Big( \hat g_{m,2} \ \hat g_{n',1} - \hat g_{m,1} \ \hat g_{n',2} \Big) \ \text{d}s_{1} \ \text{d}s_{2} \\
& + {}^{t}\hat X_{lo}^{S} \ \iint {}_{t}^{t}f_{,l} \ \hat g_{m} \ \hat g_{n} \ \Big( \hat g_{o,2} \ \ \hat g_{n',1} - \hat g_{o,1} \ \hat g_{n',2} \Big) \ \text{d}s_{1} \ \text{d}s_{2} \\
& + \int {}^{t}f \ \hat g_{m} \ \hat g_{n} \ \hat g_{n',2} \ \text{d}s_{2}
- \int {}^{t}f \ \hat g_{m} \ \hat g_{n} \ \hat g_{n',1} \ \text{d}s_{1}
\Bigg]
\end{align}
$$
Thus
$$
\begin{align}

\iint {}^{t + \Delta t}f_{i}^{S} \ \delta u_{i}^{S} \ \text{d}\ {}^{t + \Delta t}\hat S

\dot{=} \ & {}^{t + \Delta t}\lambda \ e_{ijk} \ {}^{t}\hat X_{jn}^{S} \ {}^{t}\hat X_{kn'}^{S}\ \delta \hat U_{im}^{S} \ \iint {}^{t}f \ \hat g_{m} \ \hat g_{n,1} \ \hat g_{n',2} \ \text{d}s_{1} \ \text{d}s_{2} \tag{I} \\

\\

& + \frac{1}{2} \ {}^{t + \Delta t}\lambda \ e_{ijk} \ {}^{t}\hat X_{kn'}^{S}\ \hat U_{jn}^{S} \ \delta \hat U_{im}^{S} \ 
\Bigg[
\iint {}^{t}f \ \hat g_{m} \Big( \hat g_{n,1} \ \hat g_{n',2} - \hat g_{n',1} \ \hat g_{n,2} \Big) \text{d}s_{1} \ \text{d}s_{2} \\
& + \iint {}^{t}f \ \hat g_{n} \ \Big( \hat g_{m,2} \ \hat g_{n',1} - \hat g_{m,1} \ \hat g_{n',2} \Big) \ \text{d}s_{1} \ \text{d}s_{2}
\Bigg] \tag{II} \\ 

\\

& + \frac{1}{2} \ {}^{t + \Delta t}\lambda \ e_{ijk} \ {}^{t}\hat X_{lo}^{S} {}^{t}\hat X_{kn'}^{S}\ \hat U_{jn}^{S} \ \delta \hat U_{im}^{S} \ 
\iint {}_{t}^{t}f_{,l} \ \hat g_{m} \ \hat g_{n} \ \Big( \hat g_{o,2} \ \ \hat g_{n',1} - \hat g_{o,1} \ \hat g_{n',2} \Big) \ \text{d}s_{1} \ \text{d}s_{2} \\
&+ {}^{t + \Delta t}\lambda \ e_{ijk} \ {}^{t}\hat X_{jn}^{S} \ {}^{t}\hat X_{kn'}^{S} \ \hat U_{ln}^{S} \ \delta \hat U_{im}^{S} \ \iint {}_{t}^{t}f_{,l} \ \hat g_{m} \ \hat g_{n} \ \hat g_{n,1} \ \hat g_{n',2} \ \text{d}s_{1} \ \text{d}s_{2} \tag{III} \\ 

\\

& + \frac{1}{2} \ {}^{t + \Delta t}\lambda \ e_{ijk} \ {}^{t}\hat X_{kn'}^{S}\ \hat U_{jn}^{S} \ \delta \hat U_{im}^{S} \ 
\Bigg(
\int {}^{t}f \ \hat g_{m} \ \hat g_{n} \ \hat g_{n',2} \ \text{d}s_{2}
- \int {}^{t}f \ \hat g_{m} \ \hat g_{n} \ \hat g_{n',1} \ \text{d}s_{1}
\Bigg) \tag{IV}

\end{align}
$$
[*Schwlizerhof, Karl and E. Ramm. “Displacement dependent pressure loads in nonlinear finite element analyses.” (1984). (32), P1106*]

Part IV is the boundary term of each element. The boundary terms of adjacent elements will cancel each other out during assembly. 



Suppose ${}^{t}f=0$ on the boundary of the loading area.

### Numerical Integration

#### Gauss-Legendre Quadrature

Sampling points [*Finite Element Procedures (2nd), P461-462, (5.146), (5.149)*].
$$
\forall k \in {0, 1, ..., n-1}
$$

$$
\int_{-1}^{1} \prod_{j=0}^{n-1} (r-r_{j}) \ r^{k} \ \text{d} r = 0
$$

Integration weights [*Finite Element Procedures (2nd), P456, (5.138), P462, (5.150)*].
$$
\alpha_{j} =

\frac{1}{\prod_{k=0}^{j-1}(r_{j}-r_{k}) \ \prod_{k=j+1}^{n-1}(r_{j}-r_{k})} 
\int_{-1}^{1} \prod_{k=1}^{j-1}(r-r_{k}) \ \prod_{k=j+1}^{n-1}(r-r_{k}) \ \text{d} r
$$
Calculation [*https://rosettacode.org/wiki/Numerical_integration/Gauss-Legendre_Quadrature*].

### Solver



#### Linearized Equilibrium Equation



Finite Element Procedures (2nd), P524, TABLE 6.3
$$
{}_{t}e_{ij} = \frac{1}{2} \Bigg(
\frac{\partial u_{i}}{\partial {}^{t}x_{j}} + 
\frac{\partial u_{j}}{\partial {}^{t}x_{i}}
\Bigg)
$$

#### Tangent Stress-strain Relation

Stress Vector Definition

*Finite Element Procedures (2nd), P602, (6.242)*
$$
\begin{align}
{}^{t + \Delta t}\boldsymbol{\sigma} & = \big[{}^{t + \Delta t}\sigma_{11} \quad {}^{t + \Delta t}\sigma_{22} \quad {}^{t + \Delta t}\sigma_{33} \quad {}^{t + \Delta t}\sigma_{12} \quad {}^{t + \Delta t}\sigma_{23} \quad {}^{t + \Delta t}\sigma_{31}\big]^{T} \\

& = \big[{}^{t + \Delta t}\sigma_{1} \quad {}^{t + \Delta t}\sigma_{2} \quad {}^{t + \Delta t}\sigma_{3} \quad {}^{t + \Delta t}\sigma_{4} \quad {}^{t + \Delta t}\sigma_{5} \quad {}^{t + \Delta t}\sigma_{6}\big]^{T}
\end{align}
$$
Strain Vec./include/cvm/runtime/c_runtime_api.htors Definition

*Finite Element Procedures (2nd), P602, (6.243)*
$$
\begin{align}
{}^{t + \Delta t}\boldsymbol{e} & = \big[{}^{t + \Delta t}e_{11} \quad {}^{t + \Delta t}e_{22} \quad {}^{t + \Delta t}e_{33} \quad {}^{t + \Delta t}\gamma_{12} \quad {}^{t + \Delta t}\gamma_{23} \quad {}^{t + \Delta t}\gamma_{31}\big]^{T} \\

& = \big[{}^{t + \Delta t}e_{1} \quad {}^{t + \Delta t}e_{2} \quad {}^{t + \Delta t}e_{3} \quad {}^{t + \Delta t}e_{4} \quad {}^{t + \Delta t}e_{5} \quad {}^{t + \Delta t}e_{6}\big]^{T}
\end{align}
$$


#### Interpolation

Interpolation functions.

Finite Element Procedures (2nd), P345, Figure 5.5(b)
$$
h_{1} = g_{1} - \frac{g_{9} + g_{12} + g_{17}}{2}
$$

$$
h_{2} = g_{2} - \frac{g_{9} + g_{10} + g_{18}}{2}
$$

$$
h_{3} = g_{3} - \frac{g_{10} + g_{11} + g_{19}}{2}
$$

$$
h_{4} = g_{4} - \frac{g_{11} + g_{12} + g_{20}}{2}
$$

$$
h_{5} = g_{5} - \frac{g_{13} + g_{16} + g_{17}}{2}
$$

$$
h_{6} = g_{6} - \frac{g_{13} + g_{14} + g_{18}}{2}
$$

$$
h_{7} = g_{7} - \frac{g_{14} + g_{15} + g_{19}}{2}
$$

$$
h_{8} = g_{8} - \frac{g_{15} + g_{16} + g_{20}}{2}
$$

$$
h_{j} = g_{j} \quad \forall j = 9,10,...,20
$$

$$
g_{n}(r_{0},r_{1},r_{2}) =
\begin{cases}
0 & \text{if node i is not included} \\
G_{n0}(r_{0}) G_{n1}(r_{1}) G_{n2}(r_{2}) & \text{otherwise}
\end{cases}
$$

$$
G_{ni}(r_{i}) =
\begin{cases}
\frac{1}{2}(1 + R_{ni} \ r_{i}) & R_{ni} = \pm 1 \\
1-r_{i}^{2} & R_{ni} = 0
\end{cases}
$$

Coordinate Interpolation.

Finite Element Procedures (2nd) P555 (6.127)
$$
{}^{0}x_{i} = \sum_{k=1}^{N} h_{k} {}^{0}x_{i}^{k}
$$

$$
{}^{t}x_{i} = \sum_{k=1}^{N} h_{k} {}^{t} x_{i}^{k}
$$

Isoparametric Displacement Interpolation.

Finite Element Procedures (2nd), P555, (6.128)
$$
{}^{t}u_{i} = \sum_{k=1}^{N} h_{k} {}^{t}u_{i}^{k}
$$

$$
u_{i} = \sum_{k}^{N} h_{k} u_{i}^{k}
$$




### Constitutive Relation



### Integrator



### Opensees

#### Brick Element

$$
shp =
\left( \begin{matrix}
\frac{\partial N_{1}}{\partial r_{1}}
\end{matrix} \right)
$$

