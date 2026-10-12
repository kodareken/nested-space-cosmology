# Two-room comparison on the fixed Bronnikov–Fabris background

Owning note for the 4 October 2026 external bundle stored byte-for-byte at
`lab/archive/external/nsc-two-room-comparison-2026-10-04/`.
The Swedish note, the script, and `result.json` are the original files.
`receipt.json` names their digests and the copy scope.
This page is the English record of that import. It is a scope note for
those bytes. The active magnetic programme and the papers stay as they are.

## Background

The run uses the phantom black-universe geometry already recorded for
`nsc-clock-horizon` and `nsc-dirac-tetrad` at
`7656eec34141ceca22c073bca997acefb317be05`:

\[
r=\sqrt{1+\rho^2},\qquad
A=1+3\rho+3(1+\rho^2)\bigl(\arctan\rho-\pi/2\bigr),\qquad
\beta=\sqrt{1-A}.
\]

Units are \(c=\hbar=L_{\rm throat}=1\), with benchmark mass \(m=1\).
The horizon in the saved file is \(\rho_h=1.9006916054701521\).
This is the Bronnikov–Fabris regular phantom black hole, in the black-universe
family also described by Bronnikov, Dehnen, and Melnikov. The bundle cites
arXiv:gr-qc/0511109 and arXiv:gr-qc/0611022 and does not replace their source
model. It is the imported fixed background. The magnetic Reissner–Nordström
common-action branch is a different geometry, and the mix numbers below belong
to this Bronnikov–Fabris run.

## Ray observer and angular scale

The ray observer is the unit freely falling PG normal
\(u=\partial_\tau-\beta\partial_\rho\).
For an ingoing radial null ray with conserved Killing frequency \(E\),

\[
\nu=\frac{E}{1+\beta}.
\]

The angular coefficient \(\kappa/r\) is a separate spherical scale.
In the saved stations, with \(\kappa=1\), the ray ratio
\(\nu/\nu(\rho=8)\) is \(1\) at \(\rho=8\), \(0.47282122929927417\) at the
neck \(\rho=0\), and \(0.02985652762318633\) at \(\rho=-16\).
Over the same stations the angular coefficient relative to \(\rho=8\) rises to
\(\sqrt{65}\) at the neck and then falls as the interior spheres grow.
One scale factor does not replace both maps.

## Interior Hamiltonian

On \(\rho\le 0\) the child clock is \(dT=-d\rho/\sqrt{-A}\) and
\(a_\parallel=\sqrt{-A}\).
The canonical spinor \(\chi=r\sqrt{a_\parallel}\,\psi\) removes the homogeneous
spin-connection dilution. In a constant Pauli basis the free massless block is

\[
H=\frac{k}{a_\parallel}\sigma_3+\frac{\kappa}{r}\sigma_1.
\]

Expansion anisotropy rotates that instantaneous basis. The commutator
\([H,\dot H]\) is proportional to \(k\kappa(H_\parallel-H_\perp)\sigma_2\),
so the rotation is present when \(k\kappa\neq 0\).
The saved instantaneous ratio \(\omega_{1,1}/\omega_{0,1}\) moves from
\(1.126662471472409\) at the neck to \(1.0517378926095786\) at \(\rho=-16\).

Each integrated sample starts as the positive instantaneous eigenstate at the
minimum sphere and is evolved, without reset, to \(\rho=-16\).
That initial neck mode is chosen. It is not a quantum state transported from
the parent. The saved fine negative instantaneous projections are
\(0.006208084144693428\) for \(k=1\) and \(0.008273923996503355\) for \(k=4\).
The \(k=0\) control is \(7.6512696116674235\times 10^{-34}\).

## What the run leaves open

The saved checks are finite-arithmetic, solver, and sampled-norm checks.
The file records no self-consistent state–geometry backreaction, no full
parent-to-child scattering, and no particle creation for a specified
in-vacuum and out-detector. The negative projection is a mode-basis diagnostic
on the chosen neck mode. It is not an antimatter count.

The original Swedish note has two presentation limits. Its bytes were left
as printed.

- The projection column is in percent. The refinement column beside it is the
  raw absolute difference. The raw \(k=1\) and \(k=4\) projections are the
  `result.json` values above. The printed percents \(0.62080841\) and
  \(0.8273924\) are those values times \(100\). The \(k=0\) entry
  \(7.6512696\times 10^{-32}\,\%\) is the same conversion. The refinement
  figures in that column, such as \(3.556\times 10^{-17}\) for \(k=1\), are
  not percents.
- The reported null residual is not an independent check. The ray step imposes
  the ingoing velocity \(\dot\rho=-1-\beta\), and that choice meets the radial
  null condition before the residual is sampled.

The \(k=0\) run is the analytic no-mixing case, verified by the evolution.
With the longitudinal index set to zero, \([H,\dot H]=0\): every instantaneous
Hamiltonian is proportional to the same \(\sigma_1\), so the basis does not
rotate and the prediction is no mixing. The ODE and the final projection could
still have failed that prediction. The saved value
\(7.6512696116674235\times 10^{-34}\) is a solver and basis control.

## Later use on another geometry

The same observer and the same spectral diagonalization can be applied to a
state the parent transport actually produces, on a geometry that has been
solved. The Bronnikov–Fabris projections in this bundle stay with
\(r=\sqrt{1+\rho^2}\). They are not inputs to the magnetic
Reissner–Nordström branch. Reuse on a Reissner–Nordström or coupled geometry
must recompute the observer frequency from that metric and state:
\(E/(1+\beta)\) is the static N1 benchmark on this background, and the
Bronnikov–Fabris mix numbers stay here.

## Files

| File | Bytes | SHA-256 |
|---|---:|---|
| `compare_two_rooms.py` | 14935 | `8b301708cee1818b7d3c10308abeda3d99b81c2f724fa616d0f17faeb7af4a2a` |
| `result.json` | 16895 | `8f95b61427f586e2f8dd7ddf51cfa072d7c1c13681e515a7b79f836e95388f24` |
| `forskningsanteckning.md` | 12586 | `b06d6541d7e231eefe8fdd76f9a5868f8b2ca388893bf4257ec98adb4c2f9a9b` |

`result.json` records `script_sha256` equal to the script digest above.
The script writes one new JSON and refuses to replace an existing output.
It is separate from the repository campaigns.
