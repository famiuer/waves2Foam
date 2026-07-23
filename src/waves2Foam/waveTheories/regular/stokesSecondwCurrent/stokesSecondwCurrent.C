/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "stokesSecondwCurrent.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(stokesSecondwCurrent, 0);
addToRunTimeSelectionTable(waveTheory, stokesSecondwCurrent, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


stokesSecondwCurrent::stokesSecondwCurrent
(
    const word& subDictName,
    const fvMesh& mesh_
)
:
    waveTheory(subDictName, mesh_),
    H_(readScalar(coeffDict_.lookup("height"))),
    h_(readScalar(coeffDict_.lookup("depth"))),
    omega_(readScalar(coeffDict_.lookup("omega"))),
    period_(2*PI_/omega_),
    phi_(readScalar(coeffDict_.lookup("phi"))),
    k_(vector(coeffDict_.lookup("waveNumber"))),
    K_(mag(k_)),
    Uc_(coeffDict_.lookupOrDefault<vector>("Ucurrent", vector::zero)),
    omegaAbs_(omega_ + (k_ & Uc_)),
    Tsoft_(coeffDict_.lookupOrDefault<scalar>("Tsoft", period_)),
    debug_(Switch(coeffDict_.lookup("debug")))
{
    checkWaveDirection(k_);

    // Doppler shift: absolute (fixed frame) frequency. Consistent with
    // c' = c + Uc, i.e. theta = k x - (omega + k.Uc) t + phi.
    // omegaAbs_ initialised above; report it for the user.
    Info << "stokesSecondwCurrent: intrinsic omega = " << omega_
         << ", current Uc = " << Uc_
         << ", absolute (Doppler) omega = " << omegaAbs_
         << ", apparent period = " << 2.0*PI_/omegaAbs_ << endl;

    // Warn against opposing wave-current configuration (model found to be
    // inapplicable, cf. Sun (2026), Section 6.2.2: error approx. -14 %).
    if ((k_ & Uc_) < -SMALL)
    {
        WarningIn
        (
            "stokesSecondwCurrent::stokesSecondwCurrent(...)"
        )   << endl
            << "The current opposes the direction of wave propagation "
            << "((k & Ucurrent) < 0)." << endl
            << "The stokesSecondwCurrent model severely underestimates the "
            << "wave height in this configuration and should not be used."
            << endl << endl;
    }

    if
    (
        H_/2.0 - 4.0*1.0/16.0*K_*sqr(H_)*(3.0/Foam::pow(Foam::tanh(K_*h_),3.0)
        - 1.0/Foam::tanh(K_*h_)) < 0
    )
    {
        if (debug_)
        {
            WarningIn
            (
                "label stokesSecondwCurrent::eta(point x, scalar time)"
            ) << endl << "The validity of stokes second order is violated."
            << endl << "a_1 < 4 a_2, being first and second order"
            << " amplitudes respectively." << endl << endl;
            Info << "a1 = " << H_/2.0 << " , a2 = "
                 << (1.0/16.0*K_*sqr(H_)*(3.0/Foam::pow(Foam::tanh(K_*h_),3.0)
                    - 1.0/Foam::tanh(K_*h_))) << endl;
        }
    }
}


void stokesSecondwCurrent::printCoeffs()
{
    Info << "Loading wave theory: " << typeName << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


scalar stokesSecondwCurrent::factor(const scalar& time) const
{
    scalar factor(1.0);
    if (Tsoft_ > 0.0)
    {
        factor = Foam::sin(2*PI_/(4.0*Tsoft_)*Foam::min(Tsoft_, time));
    }

    return factor;
}


scalar stokesSecondwCurrent::eta
(
    const point& x,
    const scalar& time
) const
{
    // Doppler-shifted phase: theta = omegaAbs t - k.x + phi
    scalar arg(omegaAbs_*time - (k_ & x) + phi_);

    scalar eta = (
                     H_/2.0*Foam::cos(arg) // First order contribution.
                   // Second order contribution.
                   + 1.0/16.0*K_*sqr(H_)
                   *(
                        3.0/Foam::pow(Foam::tanh(K_*h_),3.0)
                      - 1.0/Foam::tanh(K_*h_)
                    )
                   *Foam::cos(2.0*arg)
                 )*factor(time)  // Hot-starting.
                 + seaLevel_;      // Adding sea level.

    return eta;
}


scalar stokesSecondwCurrent::pExcess
(
    const point& x,
    const scalar& time
) const
{
    scalar res = 0;

    // Get arguments and local coordinate system
    scalar Z(returnZ(x));
    scalar arg(omegaAbs_*time - (k_ & x) + phi_);

    // First order contribution
    res = rhoWater_*mag(g_)*H_/2.0*Foam::cosh(K_*(Z + h_))
        /Foam::cosh(K_*h_)*Foam::cos(arg);

    // Second order contribution
    res += 1.0/8.0*rhoWater_*mag(g_)*K_*Foam::sqr(H_)/Foam::sinh(2.0*K_*h_)
        *((3.0*Foam::cosh(2.0*K_*(Z + h_))/Foam::sqr(Foam::sinh(K_*h_)) - 1.0)
        *Foam::cos(2.0*arg)
         - Foam::cosh(2*K_*(Z + h_)) + 1);

    // Apply the ramping-factor
    res *= factor(time);
    res += referencePressure();

    return res;
}


vector stokesSecondwCurrent::U
(
    const point& x,
    const scalar& time
) const
{
    scalar Z(returnZ(x));

    // Intrinsic celerity: the orbital kinematics are those of the wave in
    // the frame moving with the current, hence intrinsic omega/period are
    // used for the AMPLITUDES, while the PHASE is Doppler-shifted.
    scalar cel(omega_/K_);
    scalar arg(omegaAbs_*time - (k_ & x) + phi_);

    // First order contribution
    scalar Uhorz = PI_*H_/period_ *
                   Foam::cosh(K_*(Z + h_))/Foam::sinh(K_*h_) *
                   Foam::cos(arg);

    // Second order contribution
    Uhorz += 3.0/16.0*cel*Foam::sqr(K_*H_)*Foam::cosh(2*K_*(Z + h_))
            /Foam::pow(Foam::sinh(K_*h_),4.0)*Foam::cos(2*arg)
             - 1.0/8.0*mag(g_)*sqr(H_)/(cel*h_);

    // First order contribution
    scalar Uvert = - PI_*H_/period_ *
                   Foam::sinh(K_*(Z + h_))/Foam::sinh(K_*h_) *
                   Foam::sin(arg);

    // Second order contribution
    Uvert += - 3.0/16.0*cel*sqr(K_*H_)*Foam::sinh(2*K_*(Z + h_))
            /Foam::pow(Foam::sinh(K_*h_), 4.0)*Foam::sin(2*arg);

    // Multiply by the time stepping factor
    Uvert *= factor(time);
    Uhorz *= factor(time);

    // Generate the velocity vector: wave orbital part plus the superposed
    // background current (ramped with the same soft-start factor so that
    // wave and current develop together from a quiescent field).
    // Note "-" because of "g" working in the opposite direction
    return Uhorz*k_/K_ - Uvert*direction_ + Uc_*factor(time);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
