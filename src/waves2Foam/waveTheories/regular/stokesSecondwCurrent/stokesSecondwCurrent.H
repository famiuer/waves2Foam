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
    Foam::waveTheories::stokesSecondwCurrent

Description
    Second order Stokes theory with a superposed, depth-uniform background
    current (wave-current coupling). The wave kinematics follow the classic
    second-order Stokes solution (see stokesSecond), evaluated with a
    Doppler-shifted phase function

        theta = k x - (omega + k . Uc) t + phi

    i.e. the absolute (fixed-frame) angular frequency is

        omegaAbs = omega + (k & Uc)

    where omega is the intrinsic angular frequency of the wave (defined
    WITHOUT current, standard linear dispersion) and Uc is the background
    current velocity vector. The current velocity Uc is additionally
    superposed on the orbital velocity field. This is consistent with the
    kinematic relation c' = c + Uc (Doppler effect) used in

        Sun Haodong, "A Numerical Study on the Accuracy of Nonlinear Wave
        Generation and Wave-Current Coupling in OpenFOAM", Master Thesis,
        Kyushu University, 2026. (stokesSecondwCurrent model, Section 3.4)

    NOTE: The model is intended for waves FOLLOWING the current
    ((k & Uc) >= 0). For opposing wave-current the thesis found the model
    severely underestimates the wave height (approx. -14 %) - a warning is
    issued in that case.

    A description of the general wave framework is given in

    @article { jacobsenFuhrmanFredsoe2011,
        Author = {Jacobsen, N G and Fuhrman, D R and Freds\o{}e, J},
        title = {{A Wave Generation Toolbox for the Open-Source CFD Library:
            OpenFoam\textregistered{}}},
        Journal = {{Int. J. for Numer. Meth. Fluids}},
        Year = {2012},
        Volume = {70},
        Number = {9},
        Pages = {1073-1088},
        DOI = {{10.1002/fld.2726}},
    }

SourceFiles
    stokesSecondwCurrent.C

Author
    Derived from stokesSecond by Niels Gjoel Jacobsen (Technical University
    of Denmark). Wave-current extension following Sun Haodong (Kyushu
    University, 2026).

\*---------------------------------------------------------------------------*/

#ifndef stokesSecondwCurrent_H
#define stokesSecondwCurrent_H

#include "waveTheory.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

/*---------------------------------------------------------------------------*\
                    Class stokesSecondwCurrent Declaration
\*---------------------------------------------------------------------------*/

class stokesSecondwCurrent
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
        vector Uc_;           // Background current velocity vector
        scalar omegaAbs_;     // Doppler-shifted absolute frequency
        scalar Tsoft_;

        bool debug_;

    // Protected member functions

        scalar factor(const scalar&) const;

public:

    //- Runtime type information
    TypeName("stokesSecondwCurrent");

    // Constructors

        //- from components
        stokesSecondwCurrent
        (
            const word&,
            const fvMesh& mesh_
        );


    // Destructor

        ~stokesSecondwCurrent()
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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
