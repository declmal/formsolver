[*Finite Element Procedures (2nd), P524, TABLE 6.2*]
$$
\iiint_{{}^{0}V} {}_{0}C_{ijrs} \ {}_{0}e_{rs} \ \delta \ {}_{0}e_{ij} \ \text{d} \ {}^{0}V =
\sum_{e=1}^{E} \iiint_{{}^{0}\hat V} {}_{0}C_{ijrs} \ {}_{0}e_{rs} \ \delta \ {}_{0}e_{ij} \ \text{d} \ {}^{0}\hat V
$$

Since ${}_{0}C_{ijrs} = {}_{0}C_{ijsr}$
$$
\begin{align}

{}_{0}C_{ijrs} \ {}_{0}e_{rs}

=& \frac{1}{2} \Big( {}_{0}C_{ijrs} \ {}_{0}\hat J_{sl}^{-1} \hat h_{n,l} \hat U_{rn} +
{}_{0}C_{ijrs} \ {}_{0}\hat J_{rl}^{-1} \hat h_{n,l} \hat U_{sn} + \\
& {}_{0}C_{ijrs} \ {}_{0}\hat J_{rl}^{-1} \hat h_{n,l} \ {}_{0}\hat J_{sl'}^{-1} \hat h_{n',l'} \ {}^{t}\hat U_{kn} \ \hat U_{kn'} +
{}_{0}C_{ijrs} \ {}_{0}\hat J_{rl}^{-1} \hat h_{n,l} \ {}_{0}\hat J_{sl'}^{-1} \hat h_{n',l'} \ \hat U_{kn} \ {}^{t}\hat U_{kn'} \Big) \\

=& {}_{0}C_{ijrs} \ {}_{0}\hat J_{sl}^{-1} \hat h_{n,l} \hat U_{rn} +
{}_{0}C_{ijrs} \ {}_{0}\hat J_{rl}^{-1} \hat h_{n,l} \ {}_{0}\hat J_{sl'}^{-1} \hat h_{n',l'} \ {}^{t}\hat U_{kn} \ \hat U_{kn'} \\

=& {}_{0}C_{ijrs} \Big( {}_{0}\hat J_{sl}^{-1} \hat h_{n,l} \hat U_{rn} +
{}_{0}\hat J_{rl}^{-1} \hat h_{n,l} \ {}_{0}\hat J_{sl'}^{-1} \hat h_{n',l'} \ {}^{t}\hat U_{kn} \ \hat U_{kn'} \Big)

\end{align}
$$
Since ${}_{0}C_{ijrs} = {}_{0}C_{jirs}$
$$
\begin{align}

{}_{0}C_{ijrs} \ {}_{0}e_{rs} \ \delta \ {}_{0}e_{ij}

=&  \frac{1}{2} {}_{0}C_{ijrs} \Big( {}_{0}\hat J_{sl}^{-1} \hat h_{n,l} \hat U_{rn} +
{}_{0}\hat J_{rl}^{-1} \hat h_{n,l} \ {}_{0}\hat J_{sl'}^{-1} \hat h_{n',l'} \ {}^{t}\hat U_{kn} \ \hat U_{kn'} \Big) \\
& \Big( {}_{0}\hat J_{jl''}^{-1} \hat h_{m,l''} \delta_{ip} +
{}_{0}\hat J_{il''}^{-1} \hat h_{m,l''} \delta_{jp} + 
{}_{0}\hat J_{il''}^{-1} \hat h_{n'',l''} \ {}_{0}\hat J_{jl'''}^{-1} \hat h_{m,l'''} {}^{t}\hat U_{pn''} +
{}_{0}\hat J_{il''}^{-1} \hat h_{m,l''} \ {}_{0}\hat J_{jl'''}^{-1} \hat h_{n'',l'''} {}^{t}\hat U_{pn''} \Bigg) \delta \hat U_{pm} \\

=& \frac{1}{2} \Big( {}_{0}\hat J_{sl}^{-1} \hat h_{n,l} \hat U_{rn} +
{}_{0}\hat J_{rl}^{-1} \hat h_{n,l} \ {}_{0}\hat J_{sl'}^{-1} \hat h_{n',l'} \ {}^{t}\hat U_{kn} \ \hat U_{kn'} \Big) \\
& \Big( {}_{0}C_{ijrs} \ {}_{0}\hat J_{jl''}^{-1} \hat h_{m,l''} \delta_{ip} +
{}_{0}C_{ijrs} \ {}_{0}\hat J_{il''}^{-1} \hat h_{m,l''} \delta_{jp} + \\
& {}_{0}C_{ijrs} \ {}_{0}\hat J_{il''}^{-1} \hat h_{n'',l''} \ {}_{0}\hat J_{jl'''}^{-1} \hat h_{m,l'''} {}^{t}\hat U_{pn''} +
{}_{0}C_{ijrs} \ {}_{0}\hat J_{il''}^{-1} \hat h_{m,l''} \ {}_{0}\hat J_{jl'''}^{-1} \hat h_{n'',l'''} {}^{t}\hat U_{pn''} \Bigg) \delta \hat U_{pm} \\

=& \Big( {}_{0}\hat J_{sl}^{-1} \hat h_{n,l} \hat U_{rn} +
{}_{0}\hat J_{rl}^{-1} \hat h_{n',l} \ {}_{0}\hat J_{sl'}^{-1} \hat h_{n,l'} \ {}^{t}\hat U_{kn'} \ \hat U_{kn} \Big) \\
& \Big( {}_{0}C_{ijrs} \ {}_{0}\hat J_{jl''}^{-1} \hat h_{m,l''} \delta_{ip} + 
{}_{0}C_{ijrs} \ {}_{0}\hat J_{il''}^{-1} \hat h_{n'',l''} \ {}_{0}\hat J_{jl'''}^{-1} \hat h_{m,l'''} {}^{t}\hat U_{pn''} \Big) \delta \hat U_{pm} \\

=& \Big( {}_{0}\hat J_{sl}^{-1} \hat h_{n,l} \delta_{kr} +
{}_{0}\hat J_{rl}^{-1} \hat h_{n',l} \ {}_{0}\hat J_{sl'}^{-1} \hat h_{n,l'} \ {}^{t}\hat U_{kn'} \Big)
{}_{0}C_{ijrs} \Big( {}_{0}\hat J_{jl''}^{-1} \hat h_{m,l''} \delta_{ip} + 
{}_{0}\hat J_{il''}^{-1} \hat h_{n'',l''} \ {}_{0}\hat J_{jl'''}^{-1} \hat h_{m,l'''} {}^{t}\hat U_{pn''} \Big)
\hat U_{kn} \ \delta \hat U_{pm} 

\end{align}
$$
Thus
$$
\begin{align}
\iiint_{{}^{0}\hat V} {}_{0}C_{ijrs} \ {}_{0}e_{rs} \ \delta \ {}_{0}e_{ij}

=& \hat U_{kn} \ \delta \hat U_{pm} \iiint
\Big( {}_{0}\hat J_{sl}^{-1} \hat h_{n,l} \delta_{kr} +
{}_{0}\hat J_{rl}^{-1} \hat h_{n',l} \ {}_{0}\hat J_{sl'}^{-1} \hat h_{n,l'} \ {}^{t}\hat U_{kn'} \Big)
{}_{0}C_{ijrs} \\
&\Big( {}_{0}\hat J_{jl''}^{-1} \hat h_{m,l''} \delta_{ip} + 
{}_{0}\hat J_{il''}^{-1} \hat h_{n'',l''} \ {}_{0}\hat J_{jl'''}^{-1} \hat h_{m,l'''} {}^{t}\hat U_{pn''} \Big)
\Big|\det \boldsymbol{{}^{0}\hat J}\Big| \ \text{d}r_{1} \ \text{d}r_{2} \ \text{d}r_{3}

\end{align}
$$

Let
$$
\alpha = 3n+k
$$

$$
\beta = 3m + p
$$

$$
\phi = 2|i-j| + \bigg\lceil \frac{i+j}{2} \bigg\rceil
$$

$$
\psi = 2|r-s| + \bigg\lceil \frac{r+s}{2} \bigg\rceil
$$

Thus
$$
\begin{align}

\iiint_{{}^{0}\hat V} {}_{0}C_{ijrs} \ {}_{0}e_{rs} \ \delta \ {}_{0}e_{ij}

=& {\hat u_{a}} \ \delta \hat u_{\beta} \ \iiint

\Big( {}_{0}B^{L0}_{\psi\alpha} + {}_{0}^{t}B^{L1}_{\psi\alpha} \Big)
{}_{0}C_{\phi \psi}
\Big( {}_{0}B^{L0}_{\phi\beta} + {}_{0}^{t}B^{L1}_{\phi\beta} \Big)
\Big|\det \boldsymbol{{}^{0}\hat J}\Big| \ \text{d}r_{1} \ \text{d}r_{2} \ \text{d}r_{3}

\end{align}
$$
where
$$
{}_{0}B^{L0}_{\psi\alpha} = {}_{0}\hat J_{sl}^{-1} \hat h_{n,l} \delta_{kr}
$$

$$
{}_{0}^{t}B^{L1}_{\psi\alpha} = {}_{0}\hat J_{rl}^{-1} \hat h_{n',l} \ {}_{0}\hat J_{sl'}^{-1} \hat h_{n,l'} \ {}^{t}\hat U_{kn'}
$$

when $\psi=0, \ \alpha=3n$, which means $r=s=0$, $k=0$
$$
{}_{0}B^{L0}_{0;3n}

= \sum_{x=0}^{3} {}_{0}\hat J^{-1}_{0;x} \hat h_{n,x}
$$

$$
{}_{0}^{t}B^{L1}_{0;3n}

= 
$$




$$
{}_{0}B^{L0}_{0;3n+1} = 0
$$

$$
{}_{0}B^{L0}_{0;3n+1}
$$



