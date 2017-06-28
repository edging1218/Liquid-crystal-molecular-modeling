# Liquid Crystal Simulation
We use a Landau--de Gennes (LdG) continuum model for the Q tensor in this work. The total free energy is given as follows:
\begin{equation}
\begin{split}
f=&\int_{\mathrm{bulk}}\left(\frac{A}{2}\left(1-\frac{U}{3}\right)Q_{ij}Q_{ji}-\frac{AU}{3}Q_{ij}Q_{jk}Q_{ki}+\frac{AU}{4}\left(Q_{ij}Q_{ji}\right)^2\right)\,dV \\ 
+&\int_{\mathrm{bulk}}\left(\frac{L}{2}\frac{\partial Q_{ij}}{\partial x_k} \frac{\partial Q_{ij}}{\partial x_k} + 2q_0L\epsilon_{ikl} Q_{ij}\frac{\partial Q_{lj}}{\partial x_k}\right)\,dV \\
+&\int_{\mathrm{surf}}\left(W\left(\tilde Q_{ij} - \tilde Q_{ij}^\bot \right)^2\right)\,dS.
\end{split}
\end{equation}
