# The child glued to the parent's exterior: the neck, its opening, and what caps the mass ratio

**Corrected 7 Oct 2026 (own error): the neck-at-the-bounce number below (3.57e-23 m) is superseded. See `israel-matching-note.md`, section 0.** It used the bare density; with the geometric (torsion-included) mass the child has no effective mass at the bounce, and the mass-matched neck is a spacelike surface until t_exit/27 (R = r_g/9). The exit time, exit size and hidden-child limits below are unchanged.

7 Oct 2026. Script `junction_neck.py`, data `junction_neck.json`, figure `junction_neck.png`. Needs the restored `child_*` scripts and JSON files.
Source status: Frolov, Markov, Mukhanov, PRD 41, 383 (1990) and PLB 216, 272 (1989) were located by search but could not be read here (no open copy). Everything below is derived from the standard
Misner-Sharp mass of a closed FRW geometry; the FMM statements quoted are from memory and marked as such.

## Junction condition used

Closed FRW, areal radius `R = a sin(chi)`. The Misner-Sharp mass of the sphere at `(t, chi)` is `2 G m / c^2 = a sin^3(chi) (1 + adot^2)`, `adot = da/d(ct)`. It is geometric, so it also holds with the torsion term (which only changes the Friedmann equation).
Setting `m = M_p`, the mass seen by the parent's outside, fixes the neck: `R_j^3 = r_gp a^2 / (1 + adot^2)`. The neck is hidden inside the parent's horizon, `R_j < r_gp`, iff `a < r_gp sqrt(1 + adot^2)`.
This is the equal-mass condition. The full thin-shell (Israel) matching is in `israel-matching-note.md`.

## Results, 10 solar-mass parent (`r_gp = 29.5 km`, `G M_p / c^3 = 49 microseconds`)

- **Neck at the bounce (superseded).** The value `3.57e-23 m` came from the bare density `eps_b`. Corrected: the smallest mass-matched neck is `6.7e-23 m` at `a = 1.22 a_b`, and it is not a timelike surface (see the correction above).
- **Neck growth.** For a radiation child `R_j = 2^(2/3) r_gp^(1/3) (c t)^(2/3) = 1.59 r_g^(1/3) (ct)^(2/3)` from `t_exit/27` on, the same `2/3` power as the parent's sphere radius near the crunch, `r = 1.31 r_g^(1/3) (c tau)^(2/3)` (exact Schwarzschild interior), read backward in time. The coefficients differ by 21%, from different clocks.
  For the neck that follows the parent's thread exactly, the radius is the parent's sphere radius read in reverse (`israel-matching-note.md`, section 3).
- **Opening.** The neck reaches the parent's horizon radius at `t = G M_p / c^3 = 49 microseconds` (analytic for radiation, numerical `1.000000000000017`), at a child size `a = sqrt(r_gp r_rad) = 1.1e11 m`, which is `2.7e-7` of the child's maximum size. For an all-matter child `t = 1.22 G M_p / c^3`, `a = 1.7e11 m`.
  After that the neck is no longer hidden.
- **Staying hidden.** A child hidden through its own turnaround needs turnaround mass `<=` the parent's (`a_max < r_gp`). The largest fermion number is `3.5e58` for radiation (`2.98 N_0`) and `9.7e57` for all matter (`0.82 N_0`), with `N_0 = 1.19e58` the parent's baryons. The shear-produced `~1.3e78` does not fit.
  (From memory, FMM's semiclosed worlds stay hidden because the interior is de Sitter-like with a horizon smaller than the exterior's.)

## Reading

- For the large counts of the shear route the child is the open case: its neck opens to the parent's horizon scale after `G M_p / c^3`, and the child becomes a separate region seen from the parent's side as a small neck of mass `M_p`.
  This is the Einstein-Rosen-bridge picture of Popławski (the mother as a white hole of mass `M_p` in the child), not a hidden closed world.
- In the open case the exterior mass `M_p` does not cap the child's mass, but the neck needs a negative-energy layer (`israel-matching-note.md`, sections 2 to 4). The cap that is left is the one for the child's own collapse and recollapse, which is not computed.
- The hidden case caps the child at about the parent: no production.

## Not accounted for

Done in `israel-matching-note.md`: the neck of a closed child larger than M_p needs a negative-energy layer; what it is made of is not accounted for. The recollapse phase of the child and whether it opens a second neck. Non-radiation children for the exit time beyond the three cases listed.
