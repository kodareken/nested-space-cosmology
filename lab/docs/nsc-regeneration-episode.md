# Regeneration episode from the saved T=0.05 state

One continuation of the owned spherical Galerkin step. The source,
the feedback action, and the chart are unchanged. No restoring force,
pulse, mean subtraction, or state reset was added. The incoming gate,
a cosmological fit, and a global eternity are not requirements.

Status `FINAL`. The completed branch is two connected transfer episodes with a maintained localized structure. The joined series from T=0 is still one drift: reversal is 0, and that is not a completed renewal.

## Measured window

The coarse pilot was aimed at T=0.15. At the stored frame T=0.085 every run stopped because the both-positive expansion arc crossed into the packet. On `nf512` that arc runs from x=3.96875 to x=5.91796875. At T=0.05 its left edge was near x=5.10, on the bridge. The new left edge is 0.03125 inside the packet [0, 4), eight fine nodes, and the same coordinate on `nf256`. The arc still extends across the bridge; it has not moved wholly into the packet. The both-negative arc remains on the areal maximum, x in [1.359375, 1.62890625]. No third same-sign arc appeared. q=rQ at the areal maximum moved from 1.24054 to 1.25038 and did not turn back. r, Q, and L stayed positive. The column Gram stayed near 1e-9. The largest accepted dt*omega was 1.201, below the safety cap 1.4, so the cap did not bind and no step was rejected.

Over coordinate time 0.035 the proper clock advanced 0.10673 and the leader-window clock 0.06438. Window 0 stayed the leader. Its share moved from 0.60870 to 0.61653 and its shell content rose by 0.08117. Packet proper flux at the end is -4.98756 and the reservoir flux cancels it to about 1e-13. That cancellation is the periodic partition identity. The normal-window residual, which also carries pressure work and lapse-gradient exchange, is at most about 1e-5 of its exchange. Field energy changed by -0.28608. The coordinate Hamiltonian defect is about 2e-9 of that exchange on the confirmation step and about 6e-8 on the larger step. Four resolution pairs agree inside one percent of each claimed effect. chi_max grew from 11.148 to 50.266. This is not a global horizon.

## Handoff

Each resolution starts from the saved `dt=0.0005` final at T=0.05.
The half-step finals are siblings. They are not the initial data.

- `nf256_dtmax_0_0005` hash `f24f5dbfb5591e2943efc4d60eef0b09f573643de5ea868cea2872d46d222441`, attained T=0.085, stop `BRIDGE_ARC_ENTERED_PACKET`, steps 70, CPU 9.661763.
- `nf256_dtmax_0_00025` hash `f24f5dbfb5591e2943efc4d60eef0b09f573643de5ea868cea2872d46d222441`, attained T=0.085, stop `BRIDGE_ARC_ENTERED_PACKET`, steps 140, CPU 19.038811000000003.
- `nf512_dtmax_0_00025` hash `cd6b1f88a12635b137aaddc47444c06cc681c7568271ef989312f48f43a0649e`, attained T=0.085, stop `BRIDGE_ARC_ENTERED_PACKET`, steps 140, CPU 78.667088.
- `nf512_dtmax_0_0005` hash `cd6b1f88a12635b137aaddc47444c06cc681c7568271ef989312f48f43a0649e`, attained T=0.085, stop `BRIDGE_ARC_ENTERED_PACKET`, steps 70, CPU 38.584868.

## Protocol

Coarse budget 600 CPU seconds. Fine budget 1800 CPU seconds. Payload 64 MiB.
The pilot cap is 0.0005 and the confirmation cap is 0.00025. Each step
uses the minimum of that cap, the owned stable timestep, and the distance
to the next 0.005 frame. The safety factor 1.4 reduces the next dt.
The absolute RK4 limit is 2.8284271247461903. A step past that limit is rejected
and the last valid state is kept. Gram admissibility is 1e-8.
Declared geometric stops, evaluated on the stored frames and not on the
0.10 or 0.50 proxies: a both-positive arc entering the packet [0, 4),
a both-negative arc entering the bridge [4, 8), a third same-sign arc of
at least 4 nodes, or a turn-back of q=rQ at the areal maximum
of at least 0.001 across 2 frames.

## Forecast

```json
{
  "nf256_dtmax_0_0005": {
    "measured_step_cpu": 0.188712,
    "omega0": 939.8837108845595,
    "models": {
      "published_mean_growth": {
        "steps": 309,
        "dt_min": 1.1805e-06,
        "omega_end": 9568.150048943218,
        "growth": 23.204390026541002,
        "cpu_seconds": 67.0588092
      },
      "published_late_growth": {
        "steps": 256,
        "dt_min": 6.17515e-05,
        "omega_end": 6620.583675160021,
        "growth": 19.52182657711446,
        "cpu_seconds": 55.556812799999996
      },
      "constant_initial_dt": {
        "steps": 200,
        "cpu_seconds": 43.40375999999999,
        "growth": 0.0
      }
    },
    "pessimistic_cpu": 67.0588092
  },
  "nf256_dtmax_0_00025": {
    "measured_step_cpu": 0.188712,
    "omega0": 939.8837108845595,
    "models": {
      "published_mean_growth": {
        "steps": 432,
        "dt_min": 3.02076e-05,
        "omega_end": 9568.150048943218,
        "growth": 23.204390026541002,
        "cpu_seconds": 93.7521216
      },
      "published_late_growth": {
        "steps": 404,
        "dt_min": 0.0001015441,
        "omega_end": 6620.583675160021,
        "growth": 19.52182657711446,
        "cpu_seconds": 87.67559519999999
      },
      "constant_initial_dt": {
        "steps": 400,
        "cpu_seconds": 86.80751999999998,
        "growth": 0.0
      }
    },
    "pessimistic_cpu": 93.7521216
  },
  "nf512_dtmax_0_00025": {
    "measured_step_cpu": 0.794816,
    "omega0": 1881.8469888990094,
    "models": {
      "published_mean_growth": {
        "steps": 400,
        "dt_min": 0.00025,
        "omega_end": 4778.788300722809,
        "growth": 9.319332853294346,
        "cpu_seconds": 365.61535999999995
      },
      "published_late_growth": {
        "steps": 418,
        "dt_min": 1.04702e-05,
        "omega_end": 7562.546953174473,
        "growth": 13.909542965555442,
        "cpu_seconds": 382.06805119999996
      },
      "constant_initial_dt": {
        "steps": 400,
        "cpu_seconds": 365.61535999999995,
        "growth": 0.0
      }
    },
    "pessimistic_cpu": 382.06805119999996
  },
  "nf512_dtmax_0_0005": {
    "measured_step_cpu": 0.794816,
    "omega0": 1881.8469888990094,
    "models": {
      "published_mean_growth": {
        "steps": 243,
        "dt_min": 4.4981e-06,
        "omega_end": 4778.788300722809,
        "growth": 9.319332853294346,
        "cpu_seconds": 222.11133119999997
      },
      "published_late_growth": {
        "steps": 310,
        "dt_min": 2.18e-08,
        "omega_end": 7562.546953174473,
        "growth": 13.909542965555442,
        "cpu_seconds": 283.351904
      },
      "constant_initial_dt": {
        "steps": 200,
        "cpu_seconds": 182.80767999999998,
        "growth": 0.0
      }
    },
    "pessimistic_cpu": 283.351904
  }
}
```

## What the trajectories do

`nf256_dtmax_0_0005` runs from T=0.05 to T=0.085. Q is [0.16260481431688306, 0.2664566568533389], r is [4.216282737296435, 4.86940170486562], chi_max=50.26621726821062, proper velocity [-0.08673612115835509, 0.36772304930223104]. Leader window 0 share 0.6165300597412876, shell content 6.177244779194567. Packet flux -4.987559165620825, reservoir flux 4.9875591656208815. Proper clock 0.1067297807283187, leader clock 0.06438327627621028, field-energy change -0.2860794429811655. q at the areal maximum goes from 1.240539716809688 to 1.2503811585495024.

`nf256_dtmax_0_00025` runs from T=0.05 to T=0.085. Q is [0.16260481431506724, 0.26645665684670483], r is [4.216282737296643, 4.869401704865154], chi_max=50.266217291940706, proper velocity [-0.08673612122176236, 0.36772304900532304]. Leader window 0 share 0.6165300595501829, shell content 6.177244781693979. Packet flux -4.987559269270868, reservoir flux 4.987559269270922. Proper clock 0.10672977452848326, leader clock 0.06438327636930892, field-energy change -0.28607942497508176. q at the areal maximum goes from 1.240539716809688 to 1.2503811585494513.

`nf512_dtmax_0_00025` runs from T=0.05 to T=0.085. Q is [0.16259849778767446, 0.2664567113986368], r is [4.216282738720632, 4.8694017062229795], chi_max=50.266213742753585, proper velocity [-0.0867361744990914, 0.3677336248489724]. Leader window 0 share 0.6165300597457222, shell content 6.177244786895983. Packet flux -4.987558942658274, reservoir flux 4.987558942657826. Proper clock 0.10672977454752697, leader clock 0.06438327637953754, field-energy change -0.2860794297480709. q at the areal maximum goes from 1.2405397163987109 to 1.2503811587321414.

`nf512_dtmax_0_0005` runs from T=0.05 to T=0.085. Q is [0.16259849778930535, 0.26645671140532334], r is [4.216282738720302, 4.869401706223455], chi_max=50.266213718614665, proper velocity [-0.08673617444929466, 0.3677336250852052]. Leader window 0 share 0.6165300599367629, shell content 6.177244784395931. Packet flux -4.987558839022844, reservoir flux 4.987558839022395. Proper clock 0.10672978074736245, leader clock 0.06438327628643892, field-energy change -0.28607944774320515. q at the areal maximum goes from 1.2405397163987109 to 1.250381158732199.

The pilot decision is recorded below. Remaining cases open only when the
pilot stays in the positive chart with an admissible Gram and at least one
named physical change exceeds its predeclared floor.

```json
{
  "numerically_feasible": true,
  "named_effect_above_floor": true,
  "named_effects": [
    {
      "name": "proper_velocity_max",
      "change": 0.27956696176599005,
      "floor": 1e-08
    },
    {
      "name": "proper_velocity_min",
      "change": 0.024810513863469295,
      "floor": 1e-08
    },
    {
      "name": "proper_velocity_mean",
      "change": 0.04708287078683797,
      "floor": 1e-08
    },
    {
      "name": "observer_K_perp_max",
      "change": 0.0647229142431076,
      "floor": 1e-08
    },
    {
      "name": "observer_K_perp_min",
      "change": 0.00509897185312475,
      "floor": 1e-08
    },
    {
      "name": "observer_K_r_max",
      "change": 0.1197431804253509,
      "floor": 1e-08
    },
    {
      "name": "observer_K_r_min",
      "change": 1.1713003056630569,
      "floor": 1e-08
    },
    {
      "name": "areal_radius_min",
      "change": 0.003929802797662418,
      "floor": 1e-08
    },
    {
      "name": "areal_radius_max",
      "change": 0.00499429181391875,
      "floor": 1e-08
    },
    {
      "name": "conformal_Q_min",
      "change": 0.05691724047560967,
      "floor": 1e-08
    },
    {
      "name": "conformal_Q_max",
      "change": 0.0051732647786070785,
      "floor": 1e-08
    },
    {
      "name": "chi_min",
      "change": 0.2443072544594777,
      "floor": 1e-08
    }
  ],
  "open_remaining_cases": true,
  "stop_reason": "BRIDGE_ARC_ENTERED_PACKET",
  "comparison_target": 0.085
}
```

## Ledger

Coordinate balance compares the change in field plus gravity energy with
the integrated coordinate fieldwork. The normal-window balance compares
the change in shell content with the integral of observer boundary flux,
pressure work, and lapse-gradient exchange. A periodic flux sum of zero
is a separate partition identity and is not that balance. Mode-projector
residuals stay in the full, projected, and held-out columns.

`nf256_dtmax_0_0005` total-energy change -1.8500347920280547e-08, coordinate work -0.2860785611167556, pressure work -0.0005906149010978328, lapse exchange 0.005134452895645875. Coordinate closure within one percent: `True`. Normal-window closure within one percent: `True`.

`nf256_dtmax_0_00025` total-energy change -5.752660570124135e-10, coordinate work -0.2860792085767463, pressure work -0.0005906937146203675, lapse exchange 0.005134478654007863. Coordinate closure within one percent: `True`. Normal-window closure within one percent: `True`.

`nf512_dtmax_0_00025` total-energy change -5.749072329308547e-10, coordinate work -0.28607921335016884, pressure work -0.0005906936057535557, lapse exchange 0.005134482400096137. Coordinate closure within one percent: `True`. Normal-window closure within one percent: `True`.

`nf512_dtmax_0_0005` total-energy change -1.849213759896884e-08, coordinate work -0.286078565887805, pressure work -0.0005906147919493636, lapse exchange 0.005134456641445668. Coordinate closure within one percent: `True`. Normal-window closure within one percent: `True`.

## Comparisons

Common comparison time: `0.085`. Effect summary: `resolved`.

| Pair | Effect | Primary change | Movement | Status |
|---|---|---:|---:|---|
| time_nf256 | proper_velocity_max | 0.27956696176599005 | 2.96907998187379e-10 | resolved |
| time_nf256 | proper_velocity_min | -0.024810513863469295 | 6.340726554920906e-11 | resolved |
| time_nf256 | proper_velocity_mean | 0.04708287078683797 | 3.847703405890357e-11 | resolved |
| time_nf256 | observer_K_perp_max | 0.0647229142431076 | 6.985829970052038e-11 | resolved |
| time_nf256 | observer_K_perp_min | -0.00509897185312475 | 1.302931720847944e-11 | resolved |
| time_nf256 | observer_K_r_max | 0.1197431804253509 | 6.346090319908626e-12 | resolved |
| time_nf256 | observer_K_r_min | -1.1713003056630569 | 3.062405884435293e-10 | resolved |
| time_nf256 | areal_radius_min | 0.003929802797662418 | 2.078337502098293e-13 | resolved |
| time_nf256 | areal_radius_max | -0.00499429181391875 | 4.662936703425657e-13 | resolved |
| time_nf256 | conformal_Q_min | -0.05691724047560967 | 1.8158252679256748e-12 | resolved |
| time_nf256 | conformal_Q_max | 0.0051732647786070785 | 6.634082172496392e-12 | resolved |
| time_nf256 | chi_min | 0.2443072544594777 | 3.362016220975761e-11 | resolved |
| time_nf256 | chi_max | 39.118010861241075 | 2.3730081011308357e-08 | resolved |
| time_nf256 | weyl_C2_max | 2.313436489247834 | 2.213980110354896e-09 | resolved |
| time_nf256 | coordinate_work_integral | -0.2860785611167556 | 6.474599906769463e-07 | resolved |
| time_nf256 | proper_pressure_work_integral | -0.0005906149010978328 | 7.8813522534746e-08 | resolved |
| time_nf256 | window_normal_0 | 0.08117495107932893 | 2.499412232737086e-09 | resolved |
| time_nf256 | window_normal_1 | -0.10907077865108095 | 5.919520873476358e-09 | resolved |
| time_nf256 | window_normal_2 | 0.0033064149114795938 | 1.2652818792702192e-09 | resolved |
| time_nf256 | window_normal_3 | 0.029133172409861463 | 6.0367821852480574e-12 | resolved |
| time_nf512 | proper_velocity_max | 0.279574301183454 | 2.3623281109053096e-10 | resolved |
| time_nf512 | proper_velocity_min | -0.024779848897109 | 4.97967500567853e-11 | resolved |
| time_nf512 | proper_velocity_mean | 0.04708287079962292 | 3.847649976407297e-11 | resolved |
| time_nf512 | observer_K_perp_max | 0.06472291401876931 | 5.054756613276368e-11 | resolved |
| time_nf512 | observer_K_perp_min | -0.005091930687632775 | 1.0234143393850204e-11 | resolved |
| time_nf512 | observer_K_r_max | 0.11974715505093864 | 3.134159598516817e-13 | resolved |
| time_nf512 | observer_K_r_min | -1.1713003373047066 | 2.970024226556234e-10 | resolved |
| time_nf512 | areal_radius_min | 0.003929802012725858 | 3.2951419370874646e-13 | resolved |
| time_nf512 | areal_radius_max | -0.004994289070470614 | 4.75175454539567e-13 | resolved |
| time_nf512 | conformal_Q_min | -0.05692355708235905 | 1.6308898675987393e-12 | resolved |
| time_nf512 | conformal_Q_max | 0.005173319326815806 | 6.68654021040993e-12 | resolved |
| time_nf512 | chi_min | 0.24430582813628013 | 7.721601136267964e-14 | resolved |
| time_nf512 | chi_max | 39.11795207711774 | 2.413892019603736e-08 | resolved |
| time_nf512 | weyl_C2_max | 2.313434735265275 | 2.1893167279074532e-09 | resolved |
| time_nf512 | coordinate_work_integral | -0.286078565887805 | 6.474623638341725e-07 | resolved |
| time_nf512 | proper_pressure_work_integral | -0.0005906147919493636 | 7.88138041921474e-08 | resolved |
| time_nf512 | window_normal_0 | 0.08117494326642927 | 2.50005172119927e-09 | resolved |
| time_nf512 | window_normal_1 | -0.10907077193423298 | 5.922455859064257e-09 | resolved |
| time_nf512 | window_normal_2 | 0.0033064129363553074 | 1.268702060075455e-09 | resolved |
| time_nf512 | window_normal_3 | 0.029133179304369983 | 5.880129716473448e-12 | resolved |
| space_dtmax_0_0005 | proper_velocity_max | 0.279574301183454 | 7.339417463947395e-06 | resolved |
| space_dtmax_0_0005 | proper_velocity_min | -0.024779848897109 | 3.066496636029703e-05 | resolved |
| space_dtmax_0_0005 | proper_velocity_mean | 0.04708287079962292 | 1.2784953651312492e-11 | resolved |
| space_dtmax_0_0005 | observer_K_perp_max | 0.06472291401876931 | 2.2433828406054346e-10 | resolved |
| space_dtmax_0_0005 | observer_K_perp_min | -0.005091930687632775 | 7.041165491975324e-06 | resolved |
| space_dtmax_0_0005 | observer_K_r_max | 0.11974715505093864 | 3.974625587738201e-06 | resolved |
| space_dtmax_0_0005 | observer_K_r_min | -1.1713003373047066 | 3.164164974478467e-08 | resolved |
| space_dtmax_0_0005 | areal_radius_min | 0.003929802012725858 | 7.849365601941827e-10 | resolved |
| space_dtmax_0_0005 | areal_radius_max | -0.004994289070470614 | 2.7434481353338924e-09 | resolved |
| space_dtmax_0_0005 | conformal_Q_min | -0.05692355708235905 | 6.316606749379172e-06 | resolved |
| space_dtmax_0_0005 | conformal_Q_max | 0.005173319326815806 | 5.454820872774491e-08 | resolved |
| space_dtmax_0_0005 | chi_min | 0.24430582813628013 | 1.4263231975786272e-06 | resolved |
| space_dtmax_0_0005 | chi_max | 39.11795207711774 | 5.87841233326003e-05 | resolved |
| space_dtmax_0_0005 | weyl_C2_max | 2.313434735265275 | 1.7539825591939007e-06 | resolved |
| space_dtmax_0_0005 | coordinate_work_integral | -0.286078565887805 | 4.771049388896387e-09 | resolved |
| space_dtmax_0_0005 | proper_pressure_work_integral | -0.0005906147919493636 | 1.091484692341757e-10 | resolved |
| space_dtmax_0_0005 | window_normal_0 | 0.08117494326642927 | 7.812899660564199e-09 | resolved |
| space_dtmax_0_0005 | window_normal_1 | -0.10907077193423298 | 6.7168479667145675e-09 | resolved |
| space_dtmax_0_0005 | window_normal_2 | 0.0033064129363553074 | 1.975124286346386e-09 | resolved |
| space_dtmax_0_0005 | window_normal_3 | 0.029133179304369983 | 6.894508519650344e-09 | resolved |
| space_dtmax_0_00025 | proper_velocity_max | 0.2795743009472212 | 7.339478139134492e-06 | resolved |
| space_dtmax_0_00025 | proper_velocity_min | -0.02477984894690575 | 3.066497997081252e-05 | resolved |
| space_dtmax_0_00025 | proper_velocity_mean | 0.04708287083809942 | 1.2784419356481891e-11 | resolved |
| space_dtmax_0_00025 | observer_K_perp_max | 0.06472291396822175 | 2.0502755049278676e-10 | resolved |
| space_dtmax_0_00025 | observer_K_perp_min | -0.0050919306978669185 | 7.041168287149138e-06 | resolved |
| space_dtmax_0_00025 | observer_K_r_max | 0.11974715505125205 | 3.974619555063841e-06 | resolved |
| space_dtmax_0_00025 | observer_K_r_min | -1.171300337601709 | 3.1632411578996766e-08 | resolved |
| space_dtmax_0_00025 | areal_radius_min | 0.003929802013055372 | 7.848148797506838e-10 | resolved |
| space_dtmax_0_00025 | areal_radius_max | -0.00499428907094579 | 2.7434392535496954e-09 | resolved |
| space_dtmax_0_00025 | conformal_Q_min | -0.05692355708398994 | 6.316606564443772e-06 | resolved |
| space_dtmax_0_00025 | conformal_Q_max | 0.005173319320129266 | 5.454815626970699e-08 | resolved |
| space_dtmax_0_00025 | chi_min | 0.24430582813635734 | 1.4263567405248256e-06 | resolved |
| space_dtmax_0_00025 | chi_max | 39.11795210125666 | 5.878371449341557e-05 | resolved |
| space_dtmax_0_00025 | weyl_C2_max | 2.3134347374545916 | 1.7540072225763481e-06 | resolved |
| space_dtmax_0_00025 | coordinate_work_integral | -0.28607921335016884 | 4.773422546122674e-09 | resolved |
| space_dtmax_0_00025 | proper_pressure_work_integral | -0.0005906936057535557 | 1.0886681183278069e-10 | resolved |
| space_dtmax_0_00025 | window_normal_0 | 0.08117494576648099 | 7.812260172102015e-09 | resolved |
| space_dtmax_0_00025 | window_normal_1 | -0.10907076601177712 | 6.7197829523024666e-09 | resolved |
| space_dtmax_0_00025 | window_normal_2 | 0.0033064116676532473 | 1.9785444671516217e-09 | resolved |
| space_dtmax_0_00025 | window_normal_3 | 0.029133179310250112 | 6.8943518671815696e-09 | resolved |

## Indicator and bound

The indicator is the average of the rate's time derivative over that one realized step. It is not a bound on the metric curvature, and it is not the Euler-Lagrange substitution of p_chi_dot into R_h = chi+2. A bound would dominate the remainder on every admissible step. This record does not supply that domination. An independent metric formula can compare itself with the indicator. Agreement is a separate measurement.

The payload arrays are `frame_Q_dot`, `frame_p_chi_dot`, `frame_increment_Q_dot`,
`frame_increment_p_chi`, `frame_increment_dt`, `frame_indicator_Q_dot`, and
`frame_indicator_p_chi`, together with the coarse geometry and momenta and the
quadrature radius, conformal factor, proper velocity, and K_perp. Static lapse
and shift are rebuilt by `static_clock`. Phi columns are stored on the frames
when the 64 MiB payload allows.

## Goal

```json
{
  "goal_complete": true,
  "reasons_if_incomplete": [],
  "unresolved_comparisons": [],
  "radial_motion_alone_is_sufficient": false,
  "radial_motion_without_exchange": false,
  "structure_persists": true,
  "joined_renewal_proxy": false,
  "renewal_proxy_required": false,
  "two_transfer_episodes": true,
  "episode_one": "saved spherical feedback episode T=0 to T=0.05",
  "episode_two": "nf512_dtmax_0_00025",
  "admissible_output_input": true,
  "positive_geometry": true,
  "field_state_kept": true,
  "energy_work_within_one_percent": true,
  "meaningful_transfer": true,
  "incoming_gate_required": false,
  "lambda_cdm_required": false,
  "global_eternity_required": false
}
```

## What this does not claim

- A flux sum of zero is not a proper window balance.
- The 0.10 reversal and 0.50 share proxies are not programme requirements.
- Crossing the 1.4 safety cap is a timestep reduction, not a physical failure.
- Radial motion by itself is not regeneration and is not goal completion.
- chi+2 is not an independent metric-curvature measurement.
- The finite-step rate increment is an indicator, not a curvature bound.
- No global horizon, cosmological fit, or eternity is claimed.
- A pilot failure is a measured exit, not evidence that the regime cannot exist.

```sh
python scripts/lab.py scripts/derive_nsc_regeneration_episode.py
python scripts/lab.py scripts/derive_nsc_regeneration_episode.py --check
python scripts/lab.py -m pytest tests/test_nsc_regeneration_episode.py -q
```
