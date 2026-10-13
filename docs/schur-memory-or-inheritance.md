# Schur memory and recursive inheritance

The recursive scale map uses `unified_action.RecursiveScaleDilation`. For
spatial dimension $d$, its norm-preserving action and sheet exchange are

$$
(U_\Omega\psi)(x)=\Omega^{d/2}\psi(\Omega x),\qquad
T_\Omega=\sigma_1U_\Omega,\qquad
U_{\Omega_1}U_{\Omega_2}=U_{\Omega_1\Omega_2}.
$$

The same map scales a first-order Dirac operator, energy, cutoff and boundary
link together; the inherited spectral ratios remain fixed. Eliminating the
other sheet gives the energy-resolved relation

$$
\Gamma(E)=K(E)-\frac{B^2}{\Gamma(E/\Omega)}.
$$

At zero argument this becomes $\Gamma=K-B^2/\Gamma$. The physical branch is
continuous with $\Gamma\to K$ as $B\to0$. The executable scale map establishes
the dilation and this finite algebraic reduction. Solving the mode-resolved
tail on the warped spectrum remains a separate task. [Schur memory](schur-memory.md)
owns the intended local effective $\Lambda$ response; [inheritance](inheritance.md)
owns the regional state and energy dictionary. Original calculation records
are deposited at [Zenodo](https://doi.org/10.5281/zenodo.23090015).
