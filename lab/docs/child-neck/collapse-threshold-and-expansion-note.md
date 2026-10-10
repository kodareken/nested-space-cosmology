# What sets the collapse threshold, and why the Sun does not expand

7 Oct 2026. Computed in the session (SI units, `scipy` ODE integration); not yet a repository record.

## 1. The threshold for a neutron-star-like object

Same matter, two gravities. Matter: a cold free gas of neutrons (Oppenheimer-Volkoff 1939).

| Gravity | Result |
|---|---|
| Newtonian (weight = rest mass) | mass rises with central density and levels off at `5.73` solar masses (Chandrasekhar-type limit); a stable star exists at every density |
| General relativity (TOV) | maximum mass `0.710` solar masses at central density `4.0e18 kg/m^3`, radius `9.26 km`; denser stars are unstable |

The GR maximum is about 8 times lower. In GR the pressure that holds the star up also weighs, and
curved space strengthens the pull as the object gets compact. The push-back adds to the thing it resists.

For real neutron stars nuclear forces add push. Compactness `2GM/(Rc^2)`:
`0.34` (1.4 solar masses, 12 km) and `0.51` (2.14 solar masses, 12.4 km). The Buchdahl bound for a static
star with isotropic pressure and non-increasing density is `8/9 = 0.889`, whatever the matter.

## 2. Sun: pull against the expansion push

Round parameters: `H0 = 70 km/s/Mpc`, `Omega_Lambda = 0.7`. The outward acceleration from the
cosmological-constant term is `Omega_Lambda H0^2 r`.

| Quantity | Value |
|---|---|
| Push at 1 AU | `5.4e-25 m/s^2` |
| Sun's pull at 1 AU | `5.9e-3 m/s^2` |
| Ratio, pull over push | `1.1e22` |
| Radius where pull = push, Sun | `108 pc` |
| Same, `1e12` solar-mass galaxy | `1.08 Mpc` |
| Same, Local Group (`2e12` solar masses) | `1.36 Mpc` |
| Same, `1e15` solar-mass cluster | `10.8 Mpc` |
| Late-time expansion rate for pure dark energy | tends to `H = sqrt(Omega_Lambda) H0 = 1.9e-18 /s` (finite) |

## 3. The common requirement for any catcher

A bounce or a catch at a minimum radius needs a stress that violates the strong energy condition
(`rho + 3p < 0`): negative effective pressure. Einstein-Cartan torsion supplies it as `eps - alpha n^2`,
`p - alpha n^2`. A regular-black-hole core behaves as a cosmological constant at small radius. The
cosmological constant is the same kind of stress at large scale. The strengths differ enormously: the
dark-energy density is about `5.8e-10 J/m^3`; the stage-1 bounce state in the Popławski collapse check
is `3.9e74 J/m^3`, a ratio near `1e84`. One constant cannot do both, so a single source would need a
response that depends strongly on scale.
