# Part I: Foundations — Horizons, Thermodynamics, and Bounce Cosmology

## 1. Horizons: Causal Boundaries in Gravity

### 1.1 Schwarzschild Horizon

In general relativity the simplest black hole is the Schwarzschild black hole, characterised by mass $M$ and no spin or charge. It has a horizon at the Schwarzschild radius

$$
r_{\rm s} = \frac{2GM}{c^2},
$$

where $G$ is Newton's gravitational constant and $c$ is the speed of light [1]. When a mass $M$ is compressed inside a sphere of radius $R \le r_{\rm s}$, no light or signal can escape; the surface $r = r_{\rm s}$ becomes a one-way boundary.

A useful quantity is the average density inside the Schwarzschild radius. For a black hole of mass $M$ one finds

$$
\rho_{\rm avg} = \frac{3M}{\frac{4}{3}\pi r_{\rm s}^3} = \frac{3c^6}{32\pi G^3 M^2},
$$

showing that the average density scales like $M^{-2}$ [2]. This means very massive black holes can have extremely low average density even though the curvature near the horizon is strong.

### 1.2 Cosmological Horizons

In an expanding universe described by the Friedmann–Lemaître–Robertson–Walker (FLRW) metric, there is a scale where cosmic expansion outpaces local signals. The Hubble radius is

$$
R_{\rm H} = \frac{c}{H},
$$

where $H$ is the Hubble expansion rate. A sphere of radius $R_{\rm H}$ contains the Hubble volume; the proper radius of a Hubble sphere is $c/H$ [3].

Remarkably, the mass inside a Hubble sphere at the critical density $\rho_c = 3H^2/(8\pi G)$ satisfies $M_H \approx c^3/(2GH)$, and its Schwarzschild radius $2GM_H/c^2$ equals $c/H$. Thus the Hubble radius and the Schwarzschild radius of the Hubble mass coincide. This identity underlies black-hole cosmology models, which posit that the observable universe could be the interior of a black hole. However, such models require fine conditions and are not the standard view [4].

### 1.3 Geometry Flow and Causal Speed

A causal boundary forms when the speed of the "flow" of spacetime equals the local causal speed $c$. For a Schwarzschild black hole, the radial motion of outgoing light can be written (in the river model) as

$$
\frac{dr}{dt}_{\rm out} = c - \sqrt{\frac{2GM}{r}}.
$$

The horizon occurs where $\sqrt{2GM/r} = c$, so $dr/dt = 0$; inside this radius the geometry's inward flow exceeds the causal speed.

For the expanding universe, the proper-distance motion of incoming light satisfies

$$
\dot{D} = H D - c.
$$

The cosmic event horizon occurs when $HD = c$. Both expressions share the same structure: a horizon forms when geometric motion (infall or expansion) matches the causal speed. This observation motivates considering black holes and cosmological horizons as analogous causal membranes.

---

## 2. Black-Hole Thermodynamics

### 2.1 Bekenstein–Hawking Entropy and Hawking Temperature

Black holes obey thermodynamic laws. The Bekenstein–Hawking entropy of a black hole is proportional to the area $A$ of its horizon:

$$
S_{\rm BH} = \frac{k_B A}{4\ell_P^2},
$$

where $k_B$ is Boltzmann's constant and $\ell_P = \sqrt{\hbar G/c^3}$ is the Planck length [5]. The area law means that a black hole's entropy counts the number of microstates hidden behind the horizon; for a Schwarzschild black hole $A = 4\pi r_{\rm s}^2$, so $S_{\rm BH} \propto M^2$.

Black holes also radiate. The Hawking temperature of a Schwarzschild black hole of mass $M$ is

$$
T_{\rm H} = \frac{\hbar c^3}{8\pi G M k_B}.
$$

This temperature decreases as $M$ increases; large black holes are extremely cold.

### 2.2 Gibbons–Hawking Temperature of de Sitter Space

Spacetimes with cosmological horizons also have thermodynamic properties. In de Sitter space (a universe dominated by a positive cosmological constant $\Lambda$), an observer sees a horizon at radius $R_{\rm dS} = \sqrt{3/\Lambda}$ and experiences a thermal spectrum. The Gibbons–Hawking temperature is

$$
T_{\rm GH} = \frac{\hbar H}{2\pi k_B},
$$

with $H = \sqrt{\Lambda/3}$ [6].

### 2.3 Entropy Budgets and Information

Horizon thermodynamics suggests that black holes and cosmological horizons carry finite entropy and encode information. The Page curve and modern quantum gravity work imply that Hawking radiation can, in principle, be unitary. In other words, information may be recovered from radiation, preventing duplication across interior and exterior. This insight underlies the notion that any transition from a black-hole interior to a new universe must obey a finite entropy budget and preserve unitarity; otherwise, information would be lost or duplicated.

---

## 3. Bounce Cosmology and Maximum Compression

### 3.1 Loop Quantum Cosmology Bounce Equation

Classically, gravitational collapse leads to a singularity. In loop quantum cosmology (LQC) and related quantum gravity approaches, singularities are replaced by bounces: when the energy density $\rho$ reaches a critical value $\rho_c$, quantum effects become repulsive and the universe transitions from contraction to expansion. An effective Friedmann equation takes the form

$$
H^2 = \frac{8\pi G}{3}\rho \left(1 - \frac{\rho}{\rho_c}\right).
$$

At $\rho = \rho_c$ the right-hand side vanishes, so $H = 0$ and expansion momentarily stops. The derivative $\dot{H}$ becomes positive, leading to a bounce and subsequent expansion [7]. This equation shows that collapse halts at finite density instead of diverging to infinity.

### 3.2 Planck Density and Bounce Radius

The critical density $\rho_c$ is typically of order the Planck density, roughly $\sim 5 \times 10^{96} \; \mathrm{kg\,m^{-3}}$. If a fraction $\eta$ of a parent black hole's mass $M$ becomes the seed of a child universe and reaches this density, the bounce volume is $V_b = \eta M/(\rho_c c^2)$ and the characteristic bounce radius is

$$
r_b \sim \left( \frac{3\eta M}{4\pi \rho_c} \right)^{1/3}.
$$

For Planck-density seeds this radius is minuscule (often subatomic), necessitating subsequent inflation or rapid expansion to produce a large cosmological domain.

### 3.3 From Bounce to Expansion and Inflation

After the bounce, the scale factor $a(t)$ grows. Quantum corrections vanish and classical expansion takes over. If an inflationary phase occurs (driven by a scalar field or other mechanism), the universe can grow exponentially by a factor $e^N$ with $N \gtrsim 60$. This can stretch a tiny bounce region to cosmic scales and smooth out curvature and anisotropies.

---

## References

[1] [2] Schwarzschild radius — Wikipedia. https://en.wikipedia.org/wiki/Schwarzschild_radius

[3] Hubble volume — Wikipedia. https://en.wikipedia.org/wiki/Hubble_volume

[4] Black hole cosmology — Wikipedia. https://en.wikipedia.org/wiki/Black_hole_cosmology

[5] Black hole thermodynamics — Wikipedia. https://en.wikipedia.org/wiki/Black_hole_thermodynamics

[6] Temperature of Gravity. https://www.emergentmind.com/topics/temperature-of-gravity

[7] Quantum-Gravitational Bounce. https://www.emergentmind.com/topics/quantum-gravitational-bounce
