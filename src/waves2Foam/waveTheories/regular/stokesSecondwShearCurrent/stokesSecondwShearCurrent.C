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

#include "stokesSecondwShearCurrent.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(stokesSecondwShearCurrent, 0);
addToRunTimeSelectionTable(waveTheory, stokesSecondwShearCurrent, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


stokesSecondwShearCurrent::stokesSecondwShearCurrent
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
    Uref_(coeffDict_.lookupOrDefault<vector>("Uref", vector::zero)),
    exponent_(coeffDict_.lookupOrDefault<scalar>("exponent", 1.0/7.0)),
    dopplerModel_
    (
        coeffDict_.lookupOrDefault<word>("dopplerModel", "kirbyChen")
    ),
    UcEff_(vector::zero),
    omegaAbs_(0.0),
    Tsoft_(coeffDict_.lookupOrDefault<scalar>("Tsoft", period_)),
    debug_(Switch(coeffDict_.lookup("debug")))
{
    checkWaveDirection(k_);

    // Effective Doppler advection velocity of the sheared current
    if (dopplerModel_ == "surface")
    {
        UcEff_ = Uref_;
    }
    else if (dopplerModel_ == "depthAveraged")
    {
        UcEff_ = Uref_/(exponent_ + 1.0);
    }
    else if (dopplerModel_ == "kirbyChen")
    {
        UcEff_ = kirbyChenVelocity();
    }
    else
    {
        FatalErrorIn
        (
            "stokesSecondwShearCurrent::stokesSecondwShearCurrent(...)"
        )   << "Unknown dopplerModel: " << dopplerModel_ << nl
            << "Valid options are: kirbyChen, surface, depthAveraged"
            << exit(FatalError);
    }

    // Doppler shift consistent with c' = c + UcEff
    omegaAbs_ = omega_ + (k_ & UcEff_);

    Info << "stokesSecondwShearCurrent: intrinsic omega = " << omega_
         << ", Uref = " << Uref_
         << ", exponent = " << exponent_
         << ", dopplerModel = " << dopplerModel_
         << ", UcEff = " << UcEff_
         << ", absolute (Doppler) omega = " << omegaAbs_
         << ", apparent period = " << 2.0*PI_/omegaAbs_ << endl;

    // Warn against opposing wave-current configuration
    if ((k_ & UcEff_) < -SMALL)
    {
        WarningIn
        (
            "stokesSecondwShearCurrent::stokesSecondwShearCurrent(...)"
        )   << endl
            << "The effective current opposes the direction of wave "
            << "propagation ((k & UcEff) < 0)." << endl
            << "The underlying Doppler-shifted Stokes model severely "
            << "underestimates the wave height in this configuration "
            << "and should not be used." << endl << endl;
    }

    // Stokes second order validity
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
                "label stokesSecondwShearCurrent::eta(point x, scalar time)"
            ) << endl << "The validity of stokes second order is violated."
            << endl << "a_1 < 4 a_2, being first and second order"
            << " amplitudes respectively." << endl << endl;
        }
    }
}


void stokesSecondwShearCurrent::printCoeffs()
{
    Info << "Loading wave theory: " << typeName << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


scalar stokesSecondwShearCurrent::factor(const scalar& time) const
{
    scalar factor(1.0);
    if (Tsoft_ > 0.0)
    {
        factor = Foam::sin(2*PI_/(4.0*Tsoft_)*Foam::min(Tsoft_, time));
    }

    return factor;
}


vector stokesSecondwShearCurrent::Ucz(const scalar& Z) const
{
    // Normalised elevation above the bed: 0 at Z = -h (bed), 1 at Z = 0
    // (still water level); clipped so the current is Uref above the SWL
    // (wave crests / air side) and 0 below the bed.
    scalar s = (Z + h_)/h_;
    s = Foam::max(scalar(0), Foam::min(scalar(1), s));

    return Uref_*Foam::pow(s, exponent_);
}


vector stokesSecondwShearCurrent::kirbyChenVelocity() const
{
    // Kirby & Chen (1989):
    //   UcEff = 2k/sinh(2kh) * int_{-h}^{0} Uc(z) cosh(2k(z+h)) dz
    //
    // With sigma = z + h in [0, h] and the power-law profile
    // Uc = Uref (sigma/h)^n:
    //   UcEff = 2K * int_0^h Uref (sigma/h)^n w(sigma) dsigma
    //   w     = cosh(2K sigma)/sinh(2K h)
    //
    // evaluated with the overflow-safe identity
    //   cosh(a)/sinh(b) = exp(a - b) (1 + exp(-2a))/(1 - exp(-2b))
    //
    // (for a uniform current, n = 0, the integral returns exactly Uref).

    const label N = 1000;
    const scalar dSigma = h_/N;

    vector integral(vector::zero);

    for (label i = 0; i <= N; i++)
    {
        scalar sigma = i*dSigma;

        scalar w =
            Foam::exp(2.0*K_*(sigma - h_))
           *(1.0 + Foam::exp(-4.0*K_*sigma))
           /(1.0 - Foam::exp(-4.0*K_*h_));

        vector f = Uref_*Foam::pow(sigma/h_, exponent_)*w;

        // Trapezoidal weights
        if (i == 0 || i == N)
        {
            integral += 0.5*f*dSigma;
        }
        else
        {
            integral += f*dSigma;
        }
    }

    return 2.0*K_*integral;
}


scalar stokesSecondwShearCurrent::eta
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


scalar stokesSecondwShearCurrent::pExcess
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


vector stokesSecondwShearCurrent::U
(
    const point& x,
    const scalar& time
) const
{
    scalar Z(returnZ(x));

    // Intrinsic celerity: the orbital kinematics are those of the wave in
    // the frame advected with the (effective) current; intrinsic omega and
    // period are used for the AMPLITUDES, the PHASE is Doppler-shifted.
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

    // Wave orbital part plus the superposed sheared current, ramped with
    // the same soft-start factor so wave and current develop together.
    // Note "-" because of "g" working in the opposite direction
    return Uhorz*k_/K_ - Uvert*direction_ + Ucz(Z)*factor(time);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
