# stokesSecondwCurrent & stokesSecondwShearCurrent — wave–current coupling theories

Two new single-theory wave–current models for this fork (`of2-coupling-v0.1` line of work),
replacing the three-leg `combinedWaves` composition used in the `lc3p3wc2d`-style decks with
a self-consistent Doppler-shifted formulation.

Based on: Sun Haodong, *"A Numerical Study on the Accuracy of Nonlinear Wave Generation and
Wave-Current Coupling in OpenFOAM"*, Master Thesis, Kyushu University, 2026 (Section 3.4),
extended with the sheared-current profile of this fork's `powerLawCurrent` and the effective
advection velocity of Kirby & Chen (1989).

## Files added (untracked — carry across branch checkouts)

```
src/waves2Foam/waveTheories/regular/stokesSecondwCurrent/
    stokesSecondwCurrent.{C,H}                    # uniform current, Doppler phase
src/waves2Foam/waveTheories/regular/stokesSecondwShearCurrent/
    stokesSecondwShearCurrent.{C,H}               # power-law sheared current, Doppler phase
src/waves2FoamProcessing/preProcessing/setWaveProperties/regular/
    stokesSecondwCurrentProperties/stokesSecondwCurrentProperties.{C,H}
    stokesSecondwShearCurrentProperties/stokesSecondwShearCurrentProperties.{C,H}
add-wave-current-theories.patch                    # 2+2 lines for the two Make/files
```

## Enable & build

Your working tree was on `master` when these files were written, so the tracked `Make/files`
were NOT modified (editing them on master would block checking out `of2-coupling-v0.1`).
The new sources are untracked and survive the checkout. To enable:

```bash
cd <fork>/waves2Foam
git checkout of2-coupling-v0.1
git apply add-wave-current-theories.patch     # adds 2 lines to each Make/files
wmake libso src/waves2Foam
wmake libso src/waves2FoamProcessing          # (or ./Allwmake)
git add src/waves2Foam/waveTheories/regular/stokesSecondw* \
        src/waves2FoamProcessing/preProcessing/setWaveProperties/regular/stokesSecondw* \
        src/waves2Foam/Make/files src/waves2FoamProcessing/Make/files
git commit -m "Add stokesSecondwCurrent and stokesSecondwShearCurrent wave theories"
```

(The patch also applies on `master` with `git apply`; the sources have no dependency on
`powerLawCurrent` — the profile is re-implemented internally — so both branches compile.)

## Theory

### Common wave part
Second-order Stokes solution (identical amplitudes to `stokesSecond`, including its built-in
second-order return-flow term −gH²/(8·c·h)), evaluated with a **Doppler-shifted phase**:

```
theta    = omegaAbs·t − k·x + phi
omegaAbs = omega + (k · UcEff)          # consistent with c' = c + Uc
```

`omega`, `k` are the **intrinsic** (uncoupled) frequency/wave number from standard linear
dispersion — same convention as the thesis (wave defined without current; the Doppler effect
enters only through the phase). Orbital **amplitudes** keep the intrinsic `omega`; the
current is superposed on the velocity and ramped with the same `Tsoft` soft-start.

Note: thesis Eq. (3-40) as printed implies `omegaAbs = omega − k·Uc`, which contradicts its
own Eq. (3-29) (`c' = c + Uc`) and the dispersion relation; these implementations follow
Eq. (3-29). Verify on first run: apparent period at the inlet probe = `2π/(omega + k·UcEff)`.

### stokesSecondwCurrent
Depth-uniform current `Ucurrent` (vector). `UcEff = Ucurrent`.

### stokesSecondwShearCurrent
Power-law profile identical to this fork's `powerLawCurrent`:

```
Uc(z) = Uref · ((z − zBed)/depth)^exponent      clipped to [0, Uref]
```

(`Uc = Uref` at/above SWL, `0` at the bed, `exponent` default 1/7). The Doppler shift uses an
effective advection velocity `UcEff` selected by `dopplerModel`:

| dopplerModel | UcEff | Use when |
|---|---|---|
| `kirbyChen` (default) | `2k/sinh(2kh) ∫ Uc(z)·cosh(2k(z+h)) dz` (Kirby & Chen 1989) | general; weights the profile the way the wave actually feels it (surface-biased for short waves) |
| `surface` | `Uref` | deep-water waves, kh ≫ 1 |
| `depthAveraged` | `Uref/(exponent+1)` | long waves / shallow water |

The integral is evaluated once at construction (overflow-safe, 1000-point trapezoid);
for a uniform profile it returns exactly `Uref`. `currentSpeed()` returns `UcEff` for
GABC-type callers.

**Limitations:** irrotational-wave + shear-current superposition — wave–vorticity
interaction (shear-modified orbitals and dispersion) is neglected; valid for weak/moderate
shear. Opposing current: invalid (thesis found ≈ −14 % wave-height error); both the theory
and the Properties class warn when `(k · UcEff) < 0`.

## Usage

`constant/waveProperties.input`, sheared case (compare with the old `lc3p3wc2d` deck):

```c++
inletCoeffs
{
    waveType        stokesSecondwShearCurrent;
    depth           50.0;
    period          10.0;         // intrinsic period (uncoupled wave)
    height          6.0;
    direction       (1 0 0);
    phi             0.0;
    Uref            (0.5 0 0);    // current at the still water level
    exponent        0.142857;     // 1/7 power law (default if omitted)
    dopplerModel    kirbyChen;    // kirbyChen | surface | depthAveraged
    Tsoft           20.0;         // = 2T
    debug           yes;

    relaxationZone { ... relaxType INLET; ... }
}

outletCoeffs
{
    // Current-only target: waves absorbed, current passes through.
    // MUST mirror the inlet current, otherwise the zone kills the current.
    waveType        powerLawCurrent;
    Uref            (0.5 0 0);
    depth           50.0;
    exponent        0.142857;
    Tsoft           20.0;

    relaxationZone { ... relaxType OUTLET; ... }
}
```

For the uniform-current variant use `waveType stokesSecondwCurrent;` with
`Ucurrent (0.5 0 0);` at the inlet and `potentialCurrent` with `U (0.5 0 0);` at the outlet.

Then `waves2FoamSetWaveParameters` → `setWaveField` → solver, as usual.

## Relation to the existing combinedWaves decks

The old composition

```
combinedWaveNames (inletWave inletCurrent inletComp);
    inletWave     stokesSecond          (wave)
    inletCurrent  potentialCurrent      (current)
    inletComp     potentialCurrent      (−Q/h Stokes-transport compensation)
```

differs from the new theories in two ways:

1. **No Doppler shift** — `combinedWaves` sums the legs with the wave phase `omega·t − k·x`
   unmodified, so the imposed field's phase speed is `c`, not `c + Uc`. The interior NS
   solution then adjusts the wave riding the current away from the imposed target, which
   shows up as reflection/modulation at the inlet zone. The new theories impose the
   dispersion-consistent absolute frequency directly.
2. **Compensation leg** — `stokesSecond` (and therefore both new theories) already contains
   the second-order return-flow term `−gH²/(8·c·h)` in `Uhorz`, which is exactly the `−Q/h`
   Stokes-transport compensation. Check your old decks: an explicit `inletComp` leg on top of
   `stokesSecond` double-counts this compensation. With the new theories, **drop `inletComp`
   entirely** and do not re-add it.

## Recommended numerics (thesis campaign, uniform-current model)

Mesh on the *uncoupled* wave: 200–250 cells/λ, 20–30 cells/H — do not over-refine (the model
overestimates as the mesh refines; Medium beat VeryFine). Outlet zone ≥ 3λ (≥ 4λ for steep
waves). `Euler` time scheme (CrankNicolson diverges), `Gauss upwind` for `div(rhoPhi,U)`,
`vanLeer` for alpha, stock `fvSchemes` otherwise. Expect small positive (conservative)
wave-height bias growing with steepness (+2.2 % at ε=1/16, +4.2 % at 1/10) and current
(+7.1 % at 2.2 m/s). These numbers are for the uniform model — re-run the empty-tank
convergence study for the sheared variant before production use.
