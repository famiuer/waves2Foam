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

Class
    Foam::waveTheories::stokesSecondwShearCurrent

Description
    Second order Stokes waves superposed with a depth-sheared background
    current following a power-law profile (cf. powerLawCurrent):

        Uc(z) = Uref * ((z - zBed)/depth)^exponent,   clipped to [0, Uref]

    with Uc = Uref at (and above) the still water level and Uc = 0 at the
    bed. exponent = 1/7 gives the classic one-seventh turbulent profile.

    The wave part is the classic second-order Stokes solution (identical
    amplitudes to stokesSecond) evaluated with a Doppler-shifted phase

        theta    = omegaAbs*t - k.x + phi
        omegaAbs = omega + (k & UcEff)

    where omega is the INTRINSIC angular frequency (uncoupled wave, standard
    linear dispersion) and UcEff is the effective advection velocity of the
    sheared current felt by the wave, selected via "dopplerModel":

        kirbyChen (default):
            depth-weighted mean current of Kirby & Chen (1989),
            UcEff = 2k/sinh(2kd) * int_{-d}^{0} Uc(z) cosh(2k(z+d)) dz
            (reduces to Uref for a uniform current; weights the profile
            towards the surface for short waves)
        surface:
            UcEff = Uref (current at the still water level)
        depthAveraged:
            UcEff = Uref/(exponent + 1)

    Limitations: this is an irrotational-wave-plus-shear-current
    superposition. The modification of the orbital profiles and of the
    dispersion relation by the current vorticity (wave-shear interaction)
    is neglected; valid for weak to moderate shear. For waves opposing the
    current the underlying stokesSecondwCurrent model was found to severely
    underestimate the wave height (Sun, Kyushu University, 2026) - a warning
    is issued when (k & UcEff) < 0.

    References:
        Kirby, J.T. and Chen, T.-M. (1989). "Surface waves on vertically
        sheared flows: approximate dispersion relations". J. Geophys. Res.,
        94(C1), 1013-1027.

        Sun Haodong (2026). "A Numerical Study on the Accuracy of Nonlinear
        Wave Generation and Wave-Current Coupling in OpenFOAM". Master
        Thesis, Kyushu University. (stokesSecondwCurrent, Section 3.4)

SourceFiles
    stokesSecondwShearCurrent.C

Author
    Derived from stokesSecond (Niels Gjoel Jacobsen, DTU) and the
    powerLawCurrent theory of this fork (of2-coupling-v0.1).

\*---------------------------------------------------------------------------*/

#ifndef stokesSecondwShearCurrent_H
#define stokesSecondwShearCurrent_H

#include "waveTheory.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

/*---------------------------------------------------------------------------*\
                Class stokesSecondwShearCurrent Declaration
\*---------------------------------------------------------------------------*/

class stokesSecondwShearCurrent
:
    public waveTheory
{

protected:

    // Protected data
        scalar H_;
        scalar h_;
        scalar omega_;        // Intrinsic angular frequency (no current)
        scalar period_;       // Intrinsic period (no current)
        scalar phi_;
        vector k_;
        scalar K_;

        vector Uref_;         // Current velocity at the still water level
        scalar exponent_;     // Power-law exponent (default 1/7)
        word   dopplerModel_; // kirbyChen | surface | depthAveraged

        vector UcEff_;        // Effective Doppler advection velocity
        scalar omegaAbs_;     // Doppler-shifted absolute frequency

        scalar Tsoft_;
        bool debug_;

    // Protected member functions

        scalar factor(const scalar&) const;

        //- Local sheared current velocity at vertical position Z
        //  (Z = 0 at SWL, Z = -h at the bed; clipped power law)
        vector Ucz(const scalar& Z) const;

        //- Kirby & Chen (1989) depth-weighted mean of the current profile
        vector kirbyChenVelocity() const;

public:

    //- Runtime type information
    TypeName("stokesSecondwShearCurrent");

    // Constructors

        //- from components
        stokesSecondwShearCurrent
        (
            const word&,
            const fvMesh& mesh_
        );


    // Destructor

        ~stokesSecondwShearCurrent()
        {}


    // Member Functions
        virtual Switch conflictTSoftInitialise() const
        {
            if (Tsoft_ > 0)
            {
                return true;
            }
            else
            {
                return false;
            }
        };

        void printCoeffs();

        virtual scalar eta
        (
            const point&,
            const scalar&
        ) const;

        virtual scalar pExcess
        (
            const point&,
            const scalar&
        ) const;

        bool implementPressure() const
        {
            return true;
        };

        virtual vector U
        (
            const point&,
            const scalar&
        ) const;

        //- Effective advection speed (used by GABC-type callers)
        inline virtual vector currentSpeed()
        {
            return UcEff_;
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
