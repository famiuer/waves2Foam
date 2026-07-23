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
    Foam::stokesSecondwCurrentProperties

Description
    Pre-processing (setWaveProperties) class for the wave theory
    stokesSecondwCurrent: second order Stokes waves superposed with a
    depth-uniform background current (Doppler-shifted phase).

    The wave number and intrinsic angular frequency are computed from the
    UNCOUPLED wave parameters (period, depth) via standard linear
    dispersion, consistent with the convention in

        Sun Haodong, Kyushu University, 2026 (Master thesis, Section 3.4).

    The current is given through the vector entry "Ucurrent".

SourceFiles
    stokesSecondwCurrentProperties.C

Author
    Derived from stokesSecondProperties by Niels Gjoel Jacobsen (Technical
    University of Denmark).

\*---------------------------------------------------------------------------*/

#ifndef stokesSecondwCurrentProperties_H
#define stokesSecondwCurrentProperties_H

#include "setWaveProperties.H"
#include "stokesFirstProperties.H"

namespace Foam
{

class stokesSecondwCurrentProperties
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
        TypeName("stokesSecondwCurrentProperties");

        stokesSecondwCurrentProperties
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
