# What the child gives to the neck: energy, quanta and net baryon number swept up

7 Oct 2026. Script `neck_flux.py`, data `neck_flux.json`, figure `neck_flux.png`. Continues `child-and-neck-status-note.md` (fork 1). Needs `israel_matching.py`, `casimir_layer_check.py`, `junction_neck.py` and the restored `child_*` files.
10 solar-mass parent (`N0 = 1.19e58` baryons), N_need radiation child, one Dirac field of quanta (`g = 3.5`, mean energy `3.151 kT`).

## Method

The neck moves through the child's comoving radiation with speed `v = |a chi'|` relative to comoving observers. The exterior is vacuum, so whatever crosses the neck from the child side is absorbed by the layer. Per unit child time, with `A = 4 pi R^2`:
`dE = A (4/3) eps v gamma dt` (energy), `dN = A n_f v dt` (quanta), `dB = eta dN` (net baryon number, `eta` = net matter per quantum).
Check: along the thread path the absorbed energy from this formula agrees with the Codazzi identity of the layer (`integral of -A T_un d tau`) to `0.4%` (`0.0632` against `0.0629` at `R = 1e-3 r_g`, `0.6077` against `0.6052` at `0.1 r_g`, `1.1325` against `1.1311` at `0.5 r_g`, units `M_p`).

## Results, from the neck's start to `R = r_g`

| Neck path | energy swept (units `M_p`) | quanta swept | quanta / `N0` | half of the quanta by `R/r_g` | layer energy at the end |
|---|---|---|---|---|---|
| thread | `1.29` | `7.7e57` | `0.65` | `0.43` | `-1.91 M_p` |
| mass-matched (from `R = r_g/9`) | `3.15` | `1.75e58` | `1.47` | `0.51` | |
| layer with `P = -sigma/2`, to `0.9 r_g` | `2.35` | `1.35e58` | `1.14` | `0.57` | `-6.8 M_p` (`-6.1` at a tighter tolerance) |
| layer with `P = -sigma` | `1.7e36` | `2.2e58` | `1.8` | `0.63` | `-1.2e37 M_p` |

- The `P = -sigma/2` run leaves the thread through an unstable direction, so its end values depend on the integration tolerance by about 10% (`-6.1` to `-6.8 M_p`; child time `1.11` to `1.16 G M_p/c^3`). The earlier balanced-layers note quotes the `-6.1` run.
- Quanta swept up: `0.65` to `1.8` times the parent's baryon number, for every path. Most are taken in the last decade of radius, where the child is at `89 MeV` (thread end) or `162 MeV` (mass-matched, child time `G M_p/c^3`).
- Energy swept up for the three paths with sensible layer energies: `1.3` to `3.2 M_p`, the same order as the layer's own negative energy.

## Net baryon number and the parent's books

`B = eta N_swept`, in units of `N0`:

| `eta` | thread | mass-matched | `P = -sigma/2` layer |
|---|---|---|---|
| `6e-10` | `3.9e-10` | `8.8e-10` | `6.8e-10` |
| `9e-8` | `5.8e-8` | `1.3e-7` | `1.0e-7` |
| `1` | `0.65` | `1.5` | `1.1` |

For `eta = 6e-10` the layer takes up about `5e48` to `1e49` net baryons, `1e-9` of the parent's baryon number. The parent's books change at order one only if the child's net matter fraction is of order one. Matter and antimatter are swept up in equal numbers apart from `eta`.
In this model nothing crosses into the exterior: the exterior is vacuum Schwarzschild and the child's quanta end in the layer.

## Not accounted for

- Where the swept-up quanta and their energy go inside the layer (the layer's constituents are not fixed), and whether a layer with matter outside it would pass them on.
- The child's field content at the neck (`g = 3.5` for one Dirac field; with `g = 106.75` the quanta count rises by `(106.75/3.5)^(1/4) = 2.4` and the temperature falls).
- The stage before the neck reaches `r_g/9` on the mass-matched path (no timelike surface there).
