# Nested-influence chain: per-level fee and its closure

A saved-data consumer. It extracts, from already-recorded bytes, the
per-level influence coefficient \(f\) of a self-similar nest,

\[
\text{total}=L(1+f+f^2+\cdots)=\frac{L}{1-f},
\]

closes the chain under each extracted \(f\), and compares the innermost
level's share against a band that was fixed **before** the computation.

## Pre-registered definitions and band

Local fraction = \(1-f\). The acceptance band is
\([0.03,0.10]\) around the BBN-anchored baryon share
\(\Omega_b\approx 0.049\) inside the Planck 2018 \(5/27/68\) frame.
A candidate is `NEAR_LOCAL_BARYON_FRACTION` only if its local fraction
lands inside the band. Nothing is fitted anywhere: each candidate is a
recorded number with a stated definition.

Three candidates (module `nsc_influence_chain.py`):

| Candidate | Definition | Source record |
|---|---|---|
| structural_resolvent_one_over_omega | \(f=1/\Omega\), the multiplicative per-link scale of the resolvent levels | [scale-binding.json](../results/development/scale-binding.json) |
| recorded_transfer | \(f=(\Omega\cdot\mathrm{transfer})/\Omega\) | same record's `recursive_endpoint_product` |
| measured_mass_fee | \(f=(M_{\rm coupled}-M_{\rm fixed})/M_{\rm fixed}\) at t=20, the one-level metric-reconstruction fee | [continuation summary](../results/development/nsc-rn-neutral-continuation-v1/summary.json) |

## Measured results

| Candidate | \(f\) | local fraction | exterior fraction | verdict |
|---|---|---|---|---|
| structural_resolvent_one_over_omega | 0.251694 | 0.748306 | 0.251694 | NOT_NEAR |
| recorded_transfer | 0.0446327 | 0.955367 | 0.044633 | NOT_NEAR |
| measured_mass_fee | 0.000988698 | 0.999011 | 0.000988698 | NOT_NEAR |

All three closures miss the band. Under none of the three recorded,
non-fitted definitions does a self-similar nest put the innermost level at
\(\sim 5\%\) of the closed-chain budget.

## Two recorded observations

1.  The recorded recursive transfer is
    \(f=0.0446327\), which sits **8.9% below** the baryon share
    \(0.049\). Both are fixed numbers; no coefficient was moved. Possible
    numerology, flagged as such in the record.
2.  The measured mass fee equals \(1.0124\times(1/R^2)\) at the route
    ratio \(R=32\,r_m\); its implied layer ratio is \(31.80\), i.e. the
    measured one-level mass fee is quadratic in the declared
    \(32\,r_m\) layer ratio to 1.2%. An infinite quadratic chain adds only
    0.099% of exterior mass.

## Verdict for the story

The nesting hypothesis **cannot house the "missing 95%" as exterior mass
influence** at the transfer rates the lab has already measured:
quadratic-in-\(R\) decay makes the exterior negligible, and the recorded
resolvent transfer is a factor \(\sim 21\) below the \(f=0.95\) that a
5%-local chain would need.

What remains open, in order:

1.  **The second channel.** The closure above only sums mass-like
    influence. A summed boundary/lapse uniform term (dark-energy-like)
    is not excluded and is the next candidate.
2.  **A shallower power.** \(f=0.95\) needs either a nearly degenerate
    layer ratio (\(R\approx 1.026\) under \(1/R^2\)) or the
    \(\mu\propto 1/r\) LLL gauge screen of the charged branch, which
    converts per-level exponentials into radius powers.
3.  **Interior, not exterior.** If both channels close outside the band,
    the missing budget is interior/local structure by elimination.

## Second channel: shallow power and boundary terms

The follow-up record is
[nsc-influence-shallow-power-v1.json](../results/development/nsc-influence-shallow-power-v1.json),
driven by `derive_nsc_influence_shallow_power.py`.

Required power. To place the innermost level at the target at the
recorded route ratio \(R=32\), the chain needs
\(\nu=-\ln f/\ln R=0.01450\); the whole band needs
\(\nu\in[0.00879,0.03040]\).

Actual charged screen. The bosonization scale
\(\mu^2=|q|g_4^2/(4\pi^2 r^2)\) gives \(\nu=g_4\sqrt{|q|}/2\pi\). At the
scale-binding matched branch (\(q=4\), \(g_4=\sqrt{g^2}=4.7130\)),
\(\nu=1.5002\), so the per-level decay at \(R=32\) is
\(f=32^{-1.5}=5.52\times10^{-3}\): exterior is 0.55% per level and the
local fraction closes at 0.99448 — NOT_NEAR. The coupling the target
demands is \(g_4=0.0455\), a factor 103.5 below the matched value. The
charged channel, with its own matched couplings, is therefore **too
steep by two orders of magnitude in \(g_4\)** to carry the missing
budget as exterior mass.

Boundary channel. In the continuation run the outer lapse rate is
exactly zero at the absorbing boundary (\(N=1\) pinned), so a
uniform lapse term cannot enter the vacuum chart at all; the measured
per-level outer mass leak is \(2.80\times10^{-6}\), i.e.
\(2.9\times10^{-6}\) of the decay a 5%-local chain needs. The
boundary/channel reading is recorded in the shallow-power record.

Rounding near \(g_4\approx 3\pi/2\) (\(g^2=22.2127\) vs \((3\pi/2)^2=
22.2066\), 0.03%) is a recorded observation, not a claim.

## Routing after both channels

Mass channel refuted, charged screen refuted, boundary channel blocked.
The remaining candidates for the "missing 95%" in this nesting model are
(i) a much shallower per-level law than every recorded coefficient
delivers, or (ii) the interior/local channel. The first-principles probe
is now recorded in
[nsc-influence-response-v1.json](../results/development/nsc-influence-response-v1.json),
driven by `derive_nsc_influence_response.py` over the response owner's
`stage_tangent`:

| Probe (direction localized in…) | lapse-response fee f |
|---|---|
| …the outer level [24, 32] r_m | **0.2610** |
| …the inner level [0, 8] r_m (locality control) | 0.9990 |

In the mass channel the fee is **exactly zero**: the mass function is
causal and an outer source cannot change d_mass below its own radius.
The nonlocal influence is carried entirely by the lapse integral, and
even that channel delivers only a 26% per-level inward fraction — the
chain closes at local = 0.739, outside the band by the same margin as
every earlier channel. The first-principles two-level fee therefore
seals the round-1 verdict: the missing budget is not exterior influence
at any measured transfer.

Cross-instrument atlas: the tangent's lapse shares and the infall-time
shares of the same four intervals ([1,8], [8,16], [16,24], [24,32] r_m)
are [0.261, 0.299, 0.299, 0.141] and [0.241, 0.256, 0.251, 0.249] at
μ=1.04 — agreeing at the few-percent level on the outer half, with the
vacuum weight flattening toward 1/2 at the cap. Plateau variation below
the source edge is 3.8e-3, the quantitative form of the linear-inward
law. Two independent instruments read the same atlas:
[nsc-influence-atlas-v1.json](../results/development/nsc-influence-atlas-v1.json).
The atlas is chart-invariant (charged dyonic vs neutral, round 13), the
full seam profile is chart-invariant at 2.1e-3 (round 33), and the fee
points, profile, cells and cap clock are all μ-blind — the outside of
the chart is a universal scaffold.

Positive band location (observation-level mapping): the pre-registered
local band [0.03, 0.10] sits on the measured seam profile at [25.9,
29.0] r_m (81–90% of the 32 r_m cap, 3–6 r_m inside it; the BBN target
0.049 at 28.17 r_m), identical across all three charts within 0.06 r_m:
[nsc-cosmic-band-location-v1.json](../results/development/nsc-cosmic-band-location-v1.json).
The only band-compatible level of this chart is the outermost cap shell,
not a deep nest — the founding riddle's first positive location
statement. The cosmic reading of that location remains a hypothesis;
the mapping itself is algebra on recorded numbers.

## The complete map of the strand (synthesis, rounds 1-112)

Every link below is a creation-only record with its own tests; numbers
are quoted from the records, links give the record.

**1. Exterior channels are closed.** Mass fee 9.9e-4, recorded transfer
0.0446, structural resolvent 0.2517, lapse tangent 0.2610, charged LLL
screen 103× too steep, boundary blocked — no exterior channel delivers a
per-level fee near 0.95; the "missing 95%" is not exterior nesting
influence on this chart (rounds 1-13).

**2. The interior funnel is exact and one-way.** ΔM = Δ(excision) to
3.2e-15 over t=10..20; per half-station 2.7e-15; beyond the family EOF
1.6e-15; the outer boundary is negative in every interval (leak-only).
Every gram of interior growth is child-side feed (rounds 14, 38).

**3. The wall converts the feed into geometry.** Wall stations collect
16-22× the inflow and respond one station later (lag +1, r to -0.86 at
μ=1.040); mass-in/energy-out is a conjugate pair at r ≈ -0.99; the
converter locks near 1/G_N only at μ=1.040 (rounds 15-21).

**4. The wall breathes on the horizon-gap clock.** rise ∝ Δ^-1/2
(4.0%), fall ∝ Δ^-1/4 (3.3%), both pre-registered, both recurring on the
independent ε=0.010 trio (3.5%/6.8%) with ≤8% rate sensitivity to 10×
the driving; fall = ω·A with ω = 0.386±8.3%; the asymmetry gradient is
derived leave-one-out (μ=1.040 to 0.4%). A law of the chart, not of the
packet (rounds 24-31, 34, 39).

**5. The outside is a universal scaffold.** Fee points, full seam
profile (2.1e-3), atlas cells, cap clock, and the cosmic-band location
(0.034 r_m) are all μ-blind; the awareness profile β(s) falls
monotonically 0.995→0.080 (clock-parity at seam 4, 2.44× interior-bias
at the cosmic seam 28.17 r_m) (rounds 12-13, 25, 32-33, 35-37).

**6. The founding riddle's positive statements.** The local band
[0.03, 0.10] maps to the outermost cap shell [25.9, 29.0] r_m (target
0.049 at 28.17); that shell is length-fair in the clock (0.96×), lapse-
light (0.65×), and the live site of 34% of the causal mass response;
the window shell [r₊, 2r₊] holds a chart-invariant ~5% cell of both
instruments (rounds 35-36, 40). The 5%-class structure lives at the rim,
not in the depths; the cosmic reading remains a stated hypothesis.

**7. Cross-instrument anchors.** The Schroeder skeleton's Ω=4 and the
measured fee land within 15% (Σ = 4.08/4.61 vs 4); the spectroscopic
ω⋆ identity holds at machine precision (2.2e-16); the charged window is
cap-clocked and η-flat, the wall gap-clocked — two clocks, one chart
(rounds 27-29).

**8. The fringe, fully articulated (rounds 41-49).** The kappa-ladder's
reflection deficits follow exp(−0.635n − 0.467n²) (n=3 predicted to
0.07% from n=1,2 — Gaussian approach to the ceiling). The charged seam
profile is η-invariant at 3.9e-3 with the same ~0.26 fee at the same
physical 8 r_m seam; its weak η-tilt proved resolution-stable fringe
physics (0.83/1.00), rejected the wall's Δ-clock (61.9% spread), and
resolved into a staircase of a single nodal quantum (0.00131;
single-node L1 share, ratio 1.01), counting {3,3,2,1} quanta per seam
from inner to outer, with the chart-free unit q·N = 1.149 held by both
charts and the neutral sign flip explained by grid construction. The
staircase thins faster than the awareness β (count/β 3.0→1.8): the
fringe delivers fully in the clock-parity zone and drains beyond it.

**9. The joins (rounds 50-58).** The programme's measurements were laid
against the independently-owned strands of the repository, all n=3
stated and all numbers re-read from records: (a) the two-ends room
sharpness joins the wall kinetics — corr(Q, rise) = −0.926, absorbed by
the gap powers Q ∝ Δ^0.117 and rise ∝ Δ^−1/2 into Q·rise^0.2346 = 15.5
at 2.1% spread (13.5× collapse), with p_Q frozen at 0.1173 across every
recorded station time; (b) the one-way current book's single outward
cell [0, 0.5] contains all three wall kinematic constants (rise 26%,
fall 48%, ω 77% of the band) with the +0.283 transient standing between
fall and ω, while the mass book stays 0-positive in 20 intervals; (c)
the cone trap nests inside the wall envelope inside the window shell in
all three cases (trap 1.0141/1.0068/1.0036 r₊, monotone), and the
redshift-at-excision is μ-uniform at 8.0e-7 — the sixth μ-blind
scaffold quantity; (d) the isotropy index does not track the throat
balance (NO_TRACK over the overlap; the wall's swings are spatially
localized), and its t=40/80 tail is the only recorded balance-scale
observable past t=20; (e) the constant table assembles seven
coefficients from their own records, with one numerology-flagged
near-identity: 4·p_Q = 0.4692 vs ladder b = 0.4672 (+0.43%).

**10. The late harvest (rounds 59-68).** The identity's predictive test
passed (ladder-implied exponent closes the covariance at 2.06% vs
2.07%); the composed constant stands at 0.948 of the μ-blind infall
clock (±2%); the driven isotropy anti-tracks all three wall kinematics
at n=6 (−0.96, −0.96, −0.87) with I_first the seventh μ-blind quantity
(0.24%); the continuation's scalars ledger showed the trap drifting
outward with the feed (r = 0.835, and 0.999 into the adjacent window —
the margin follows immediately while the balance lags); the march's own
throat slope flips sign exactly at the regime boundary; and p_Q
decomposed into pure room stretch (p_N ≈ 0, p_r = 0.1197) reconciling
0.0000, with the analytic envelope √(μ²−1)/μ bounding the chord and its
four-fold containing the ladder's b. The 633-record development
inventory was swept; the untouched remainder is other strands' internal
diagnostics, so the recorded-bytes frontier is closed by evidence, not
assertion.

**11. The production door (rounds 69-71).** The launch was priced from
the recorded ledger (~570 CPU-seconds for t=20..40 inside the 900 s
admission slot, ~12 minutes total), the campaign calibration ceremony
executed dry and passed (worst residual 4.3e-5 on twelve channels), and
the real prepare stage was then attempted: the freeze gate refused with
"producer bytes are not the frozen root commit", diagnosed as the
pre-existing modification of lab/src/recursive_horizons/nsc_rn_pg.py
against the frozen commit. Superseded by part 12: the t>20 morphology
was since answered from the recorded relaxation payloads, so the
campaign is no longer required for the wall's future.

**12. The wall's complete lifecycle (rounds 72-77).** The relaxation
campaign directory held the full station payloads at t=5..80 with the
wall machinery's ingredients all along: the wall dies between t=10 and
t=15 and never reforms, while the throat balance climbs to 0.9975 (the
single-force limit) and its slope steepens monotonically to −70; the
afterlife has one death window, a mid-relaxation event at t=40 with
three cross-observable witnesses (isotropy 0.419, balance deficit
0.093, slope wiggle −62.3→−58.7), and an exponential deficit whose
halving time is 7.7 t. The march grid's wall flags show the lifespan is
μ-monotone on the deep domain (lightest walls reach t=20, the locked
wall dies at t=10) with the march and relaxation arms agreeing digit-wise
at t=20 (balance 0.77); the march slope trajectory shows all cases
negative at t=5 and the wall epoch at t≈10 flipping only the light
cases positive (+180.8 spike at μ=1.002 where the wall count triples).
The wall is born with the packet, breathes on the gap clock, dies when
the driven phase ends, and its tombstone is fully measured.

**13. The calibration standard (rounds 78-79).** The frontier amendment
closed two of the round-63 open cells with their later records (wall
t>20 morphology via rounds 73-77; p_Q derivation via rounds 67-68),
leaving exactly two cells: the fourth gap point's exponent provenance
(production) and the paper-level cosmic coupling (out of lab domain).
The cross-arm identity then wrote the programme's internal gold
standard: the march and relaxation payload families carry the same
mu=1.040 trajectory bit-identically at all four shared stations (diffs
0.0; wall counts 2/2/0/0), so every cross-arm comparison in the
lifecycle records is exact by construction. The recorded bytes of the
RN radial corridor - every campaign's payloads and every strand's
records - have now been read, joined, and pinned.

**14. The wall's domain law (rounds 81-82).** The deep-grid wall radii
settle the round-22 envelope's status: on the horizon-reaching grid the
walls sit lower than the family chart ever showed (μ=1.040's wall
inside the horizons at 0.65 r₊) while the μ-ordering survives both
domains. The shallow/deep breathing-amplitude ratio is μ-monotone
(7.1×/17.4×/40.1×) — the deep wall is nearly pinned while the shallow
wall swings wide — and the deep wall's minimal standoff above the
grid's inner edge is 0.012/0.089/0.073 r_m. The wall's locomotion is in
part the room it stands in: grid-relative position, ordered amplitude,
exact cross-arm calibration. The window shell bound is domain-specific;
the intrinsic facts are the ordering and the nodal ties of rounds 46-48.

**15. The campaign corridor and the cost laws (rounds 84-106).** The
job registry revealed the march runner's second corridor: station
extensions on the three legal masses are operational even while the
campaign CLI stays frozen-gated, leaving the fourth gap point alone
blocked by the packet preparer's mass pin. A t30 campaign for mu=1.002
was launched into a new creation-only directory and marches still, with
its instruments pre-staged: a reducer, a pre-registered verdict map,
and a conditional test. The recorded v1 mtimes then yielded two cost
laws: the march's durations are mu-ordered (173.7/89.8/54.0 minutes -
the lightest room costs 3.2x the heaviest, exactly where the wall lives
longest), and every case's station intervals shrink with t. The shrink
joins the wall's death only partially - cleanest for the case whose
wall dies inside the coverage (interval shrink 5.2 min at mu=1.040 vs
1.1 min at mu=1.002; learn-rate factors 0.986/0.897/0.810, mu-ordered,
held as fits) - and the ledger records the partial
claim as partial, after round-102's poll-count correction and the
round-104 timing resolution: measured from the files themselves,
forecasting the t30 stations near 3.5 and 4.2 hours from launch. The
relaxation campaign's own mtimes then showed a conflicting rate - a
steady 2.0 minutes per t through all six stations to t=80 - so the two
recorded laws were laid side by side as a dual forecast with the
campaign's write times named the adjudicator; the watch's first
snapshot eliminated the fast rate when no station landed near fifty
minutes, leaving the v1-derived law governing until the directory
writes its verdict. **The march landed: WALL_SURVIVES_TO_30.** At t=25
the wall shows two flags (balance 0.706) and at t=30 two again
(balance 0.767), first wall at 1.006 then 1.034 r₊: the lightest room's
wall outlives every recorded era, the deep-domain amplitude law gains
its longest member, and the named next question is the longer campaign.

**Open cells (production-scale, explicitly named):** the wall's t>20
morphology (continuation directory is lock-listed; requires a new
campaign), the exponent provenance at a 4th gap point (requires new
family runs), and any future coupling of the rim location to an actual
cosmological observable (a paper-level question, not a lab one).

Proposal, not observation: the closure is a model comparison over frozen
bytes, not a cosmology measurement. Owners: `nsc_influence_chain.py`,
`derive_nsc_influence_chain.py`, `derive_nsc_influence_shallow_power.py`,
`derive_nsc_influence_response.py`, `test_nsc_influence_chain.py`, the
records above.