# FGC-1-HYP1-DOM2-MODE1: pointwise REF1 branch-mode controls

MODE1 is a pointwise exact rational control artifact.  It compares the direct
REF1 branch principal matrix against the MHG1 standard first-order principal
matrix on the FO1 flat control, the FO1 activated off-shell nonzero-acceleration
control, and the solved COMP1 exact-root two-jet.  It serializes exact RREF bases at
the tilde and hat auxiliary roots.  At the COMP1 datum it also derives the
exact fixed-jet acceleration root from the complete REF1 acceleration block,
substitutes it into all six REF1 residual rows, and checks strict inclusion in
the inherited QIFT1 acceleration box.  This is not a nonlinear acceleration
map, an initial slice, or an evolution result.

The artifact does not prove a uniform-domain eigenframe or symmetrizer, any
constraint propagation, an EFT envelope, an IBVP, evolution, collapse, or
defocusing.  Optional polynomial/quotient-atlas/zero-speed controls are only
recorded when the core exposes the corresponding exact identity.

It additionally records a universal formal REF1 connection-variation/MHG
projector principal identity and pointwise physical metric/regulator quotient
charts at the solved COMP1 root.  These are not uniform-domain statements.

Reproduce with:

```bash
python3 scripts/reproduce_fgc_hyp1_dom2_mode1.py \
  --config configs/fgc/fgc-1-hyp1-dom2-mode1.toml \
  --output results/fgc-1-hyp1-dom2-mode1.json
```
