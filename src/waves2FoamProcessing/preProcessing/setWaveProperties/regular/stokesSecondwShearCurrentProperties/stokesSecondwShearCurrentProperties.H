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
    Foam::stokesSecondwShearCurrentProperties

Description
    Pre-processing (setWaveProperties) class for the wave theory
    stokesSecondwShearCurrent: second order Stokes waves superposed with a
    depth-sheared power-law background current (Doppler-shifted phase with
    a Kirby & Chen (1989) effective advection velocity by default).

    The wave number and intrinsic angular frequency are computed from the
    UNCOUPLED wave parameters (period, depth) via standard linear
    dispersion. The current profile is given through "Uref" (velocity at
    the still water level), "exponent" (default 1/7) and the optional
    "dopplerModel" (kirbyChen | surface | depthAveraged).

SourceFiles
    stokesSecondwShearCurrentProperties.C

Author
    Derived from stokesSecondProperties (Niels Gjoel Jacobsen, DTU).

\*---------------------------------------------------------------------------*/

#ifndef stokesSecondwShearCurrentProperties_H
#define stokesSecondwShearCurrentProperties_H

#include "setWaveProperties.H"
#include "stokesFirstProperties.H"

namespace Foam
{

class stokesSecondwShearCurrentProperties
:
    public setWaveProperties
{
private:

    scalar period_;
    scalar omega_;
    scalar depth_;

    stokesFirstProperties sfp_;

public:

    //- Runtime type information
        TypeName("stokesSecondwShearCurrentProperties");

        stokesSecondwShearCurrentProperties
        (
            const Time&,
            dictionary&,
            vector,
            bool
        );

        // Method

        virtual  void set( Ostream&);

};

}

#endif

// ************************************************************************* //
