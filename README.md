# Nested-Space Cosmology

The whole proposal is one recursive spectral equation, not a list of separately fitted laws:

$$
\boxed{
Z_{\Lambda}^{\mathrm{one}}(\mathbb D_\Theta)
=
\int \mathcal D\varphi\,
Z_\Lambda\!\left(
e^{-\varphi/2}\mathbb D_\Theta e^{-\varphi/2}
\right),
\qquad
\mathcal T_\Theta^*\mathbb D_\Theta=\mathbb D_\Theta.
}
$$

Its present concrete regulated action uses the heat-spectral profile
\(\operatorname{Tr}e^{-\mathbb D_\Theta^2/\Lambda^2}\) plus the fermionic
quadratic form. Every named local quantity is supposed to be a projection of
this one operator: constants, particle poles, wave residues, nuclear energy,
unresolved (“dark”) metric response, and the black-to-child boundary. The
configuration of a room can change. The law, the dimensionless coefficient
vector \(\Theta\), and the boundary kernel must not.

This repository is a conjectural research programme. It does **not** prove Nested-Space Cosmology, a final inheritance scale \(\zeta\), a particle spectrum, a dark-sector fit, or that our universe is inside a black hole.

The conceptual thesis, in one chain:

```text
spaces form inside spaces
  -> local constants are formation parameters of one room
  -> singularity and infinity mark failed local charts, not completed objects
  -> parent collapse and child expansion are two sides of one finite transition
  -> unresolved rooms act locally through the same kernel
  -> the next generation repeats the same law
```

## Authorship and AI contribution

**Douglas Ek** is the accountable author: research direction, conceptual synthesis, scientific decisions, and publication responsibility.

The public documents, derivations, scripts, and compact results were developed in iterative collaboration with **ChatGPT** and **OpenAI Codex**. The AI collaborator drafted formal language, code, literature comparisons, and adversarial checks. Every AI contribution is candidate work until it matches explicit mathematics, committed computation, cited sources, and the JSON nonclaims. AI credit records material assistance; it does not transfer scientific accountability.

## Current frontier

**The current frontier is the unsolved energy-resolved recursive child tail.**

Anomaly compensation removed the determinant-only scale roots. Adding an independently weighted geometric action is forbidden double-counting. The missing owner is the child’s own next-child response. At a fixed common dimensional energy, and with \(b=B_{\mathrm{dim}}/\Lambda_{\mathrm{parent}}\), the parent block is

$$
\Gamma_p(x)
=
K_p(x)
-
\frac{1}{\Omega}\,
b\,
\Gamma_c(x/\Omega)^{-1}
b^{\dagger}.
$$

The unweighted display \(\Gamma(x)=K(x)-B^{\dagger}\Gamma(x/\Omega)^{-1}B\) is the same equation after a symmetric link convention; the two normalizations cannot silently share a numerical \(B\). \(\Omega\) is not selected. Dilation alone does not derive \(\zeta=\Omega^2\). The mode-resolved tail is not solved on the warped transition spectrum. Physical \(\zeta\) is not derived. See [docs/current-result.md](docs/current-result.md).

## What is imported

These are launch surfaces, not project novelty. Primary sources are collected
in [docs/prior-art-and-open-claim.md](docs/prior-art-and-open-claim.md) and the
technical paper's references.

| Subject | Status in this repository |
|---|---|
| Single-electron interference, QFT field excitations, Bell nonlocality, PBR | Experimental and theoretical constraints on any particle ontology |
| Skyrme/BPS solitons and nuclear binding | Construction pattern for stationary charge sectors |
| Einstein–BPS compact-star branches | Gravity benchmark from nuclear-fixed parameters |
| Exact regular black-universe geometry | Limiting trapped-to-expanding four-geometry |
| RS2 projected Weyl term as outside response | Structural realization of unresolved-room geometry |
| Two-wall shadow matter and the \(3/4\) lensing discriminator | Other-room gravity acting on the local sheet |
| Chamseddine–Connes spectral action | Local-room Standard Model plus gravity scaffold |
| Local Higgs–dilaton scale anomaly | Common local Weyl mode, not two independent sheets |
| Doubled spectral FLRW / relative-metric interaction | Leading two-geometry interaction |
| Nonlocal spectral propagation wall | High-resolution cutoff of local boson propagation |
| Wetterich exact effective-action flow | Resolution-flow mathematics |
| Baby-universe / cosmological-natural-selection literature | Prior hypothesis class, not an observation |

## What this repository derives

The repository binds imported components to one operator, derives identities that those components do not already state, and records numerical diagnostics with explicit nonclaims.

Derived here, always conditional on the stated premises:

- one immutable \(\Theta\) with no hard-coded sector masses;
- nested \(L0\)–\(L6\) self-equality reducing Skyrme scale coefficients to one overall normalization;
- the dictionary that *names* local constants as Hessian/spectral outputs of \(S_{\mathrm{one}}\);
- the two-sheet representation of a noninjective areal-radius black-universe chart;
- the transition-fixed Gauss–Bonnet coupling \(\alpha_{\mathrm{GB}}/L_\star^2=1015/144\) and the on-room vacuum-form component \(108/1015\);
- the scalar recursive outside kernel \(\Gamma=K-b^2/\Gamma\);
- the exponential heat profile from composition of resolution steps, as a heat-kernel diagnostic rather than a complete anomaly;
- the local identity that one off-diagonal \(\Phi\) is mass gap, Schur self-energy, and relative-metric link;
- child-side orientation of the source sign, later corrected by the exact child asymptotic;
- a global horizon-penetrating spectral foliation; throat maps \(N_{\mathrm{parent}}(E)\), \(N_{\mathrm{child}}(E)\), and \(\Phi_{\mathrm{throat}}=R(E)\) remain unevaluated declarations;
- the unitary dilation that makes \(\mathcal{T}_\Theta\) executable; \(\zeta=\Omega^2\) is not derived by that dilation alone;
- the unit-consistent additive-gap argument \(q_j=\lambda_j/\zeta+\mu^2\) and a sign theorem for that finite, fixed-geometry family.

Numerical diagnostics that **must not** be read as a final \(\zeta\) or a particle-to-cosmology number are listed in [docs/current-result.md](docs/current-result.md).

## What the anomaly invalidated

Compensating the pure cutoff anomaly of the invariant partition function **invalidated the determinant-only inheritance-scale roots**. Those historical numbers also used the inconsistent mass conversion \((\lambda_j+\mu^2)/\zeta\). The declared additive-gap operator requires \(q_j=\lambda_j/\zeta+\mu^2\). After that unit correction, the same inherited cutoff subtraction has a strictly negative derivative for every \(\zeta>3\pi/2\) in the finite, fixed-geometry family. Those earlier numbers remain useful diagnostics. They are not stationary points of the one equation. Raw finite-determinant invariance does not set the regulated relative Hessian to zero.

Related routes already rejected on their own JSON nonclaims, and not repaired by retuning \(\Phi\) or adding a second geometric weight:

- \(\Lambda L_\star=1\) as the physical child-matched gap (the exact child vacuum requires \(\zeta>3\pi/2\));
- the lowest-mode heat-only scale minimum (the angular tower removes the root);
- independently weighted Einstein–Gauss–Bonnet or heat geometry as the missing scale owner (double-counting);
- pure critical five-dimensional Einstein–Gauss–Bonnet as the complete interface law (time-principal rank \(5/10\); no scanned diagonal normalization restores real characteristics);
- the standalone bosonic heat-trace graviton covariance as a positive Stieltjes propagator;
- the inverse relative heat Hessian as Osterwalder–Schrader positive metric covariance;
- constant two-wall shadow matter as the dominant kiloparsec gravitating mass.

## Map

```mermaid
flowchart TB
  subgraph imported["Imported launch surfaces"]
    BU["Black-universe geometry"]
    SA["Spectral action and scale anomaly"]
    SK["Skyrme / BPS nuclear pattern"]
    SH["Two-wall shadow-matter discriminator"]
  end

  subgraph derived["Repository derivations"]
    OP["One recursive spectral operator"]
    TS["Two-sheet carrier and GB identities"]
    RK["Recursive outside kernel"]
    ID["Local Phi = gap = outside self-energy"]
    DL["Unitary dilation T_Theta"]
  end

  subgraph invalidated["Invalidated as physical owners"]
    Z1["Direct one-scale gap 1006/1015"]
    HT["Heat-only angular zeta"]
    DT["Determinant-only zeta after anomaly compensation"]
    UC["Unit-corrected finite-family determinant"]
    EG["Diagonal critical EGB interface"]
    SW["Constant two-wall shadow dominance"]
  end

  subgraph frontier["Current frontier"]
    TL["Unsolved energy-resolved recursive child tail"]
    OZ["Open: physical zeta, spectrum, cosmology"]
  end

  imported --> OP
  OP --> TS
  OP --> RK
  RK --> ID
  ID --> DL
  TS --> invalidated
  DL --> TL
  DT --> UC
  UC --> TL
  TL --> OZ
```

## Reproduce the public result

The complete current dependency closure is small enough to recompute on a
laptop. The v0.1.0 historical 58-record chain remains byte-for-byte. A 59th
follow-up, the unit-closure check, is reproduced by comparing every exact and
numeric field. The reference historical run regenerated those 58 artifacts
byte-for-byte in about four minutes:

```bash
python3 -m pip install -e ".[paper]"
make reproduce
make paper-check
```

Run `make verify` for the integrated manifest, focused-test, reproduction, and
paper gate. Regeneration happens in an isolated temporary workspace and never
rewrites the tracked result files. See
[docs/reproducing.md](docs/reproducing.md) for exact and cross-platform modes.

## Documents

| File | Role |
|---|---|
| [THEORY.md](THEORY.md) | Public theory with labeled evidence status |
| [AGENTS.md](AGENTS.md) | Short working orientation for the one-equation thesis |
| [docs/prior-art-and-open-claim.md](docs/prior-art-and-open-claim.md) | Imported results and the one open joint claim |
| [docs/current-result.md](docs/current-result.md) | Present compact checkpoint and nonclaims |
| [docs/nsc-closure-verification-2026-09-07.md](docs/nsc-closure-verification-2026-09-07.md) | Unit-closure follow-up and remaining equations |
| [docs/reproducing.md](docs/reproducing.md) | How to inspect and regenerate compact results |
| [results/README.md](results/README.md) | Category map for the 58-record chain plus the 59th follow-up |
| [CONTRIBUTING.md](CONTRIBUTING.md) | Claim discipline and contribution rules |
| [paper/nested-space-cosmology.md](paper/nested-space-cosmology.md) | Technical paper source |
| [paper/nested-space-cosmology.pdf](paper/nested-space-cosmology.pdf) | Reproducibly generated working preprint |

## How to read a claim

Statements in this repository are labeled in [THEORY.md](THEORY.md):

- **Postulate** — a model requirement, not a theorem.
- **Imported result** — established literature or observation, used as a launch surface.
- **Repository derivation** — an identity or closure obtained here from stated premises.
- **Numerical diagnostic** — a computed number whose JSON record forbids promotion.
- **Open prediction** — a quantity that would count as a result if derived without using the compared datum.
- **Interpretive hypothesis** — a picture that is not itself an action, spectrum, or observation.

If a compact JSON `nonclaims` flag is `false`, public prose must not say the corresponding statement has been shown.

## Citation and licence

Citation metadata are in [CITATION.cff](CITATION.cff). Douglas Ek is the
accountable author; the AI contribution is disclosed above and in the paper.
Code is MIT-licensed, original prose and figures are CC BY 4.0, and cited
third-party work retains its original terms. See [LICENSE](LICENSE).

Security reports use [SECURITY.md](SECURITY.md); scientific disagreements and
reproduction failures belong in the ordinary issue tracker.
