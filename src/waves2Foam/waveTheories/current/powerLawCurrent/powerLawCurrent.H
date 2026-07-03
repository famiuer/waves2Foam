/*---------------------------------------------------------------------------*\
    powerLawCurrent — depth-sheared current with a power-law profile:

        U(z) = Uref * ((z - zBottom)/depth)^exponent ,  zBottom = localSeaLevel - depth

    clipped to [0, 1] in the normalised elevation, i.e. U = Uref at the still
    water level and above, U = 0 at the bed. exponent = 1/7 gives the classic
    one-seventh power law for a turbulent open-channel / tidal current.

    Added for the OC4 LC3.3 wave+current cases (lc3p3wc2d). Mirrors
    potentialCurrent in every other respect (eta = localSeaLevel, hydrostatic
    excess pressure, Tsoft sine ramp).

    NOTE: assumes the mesh vertical axis is global z (as all cases in this
    repository; waves2Foam itself is used with z vertical here).
\*---------------------------------------------------------------------------*/

#ifndef powerLawCurrent_H
#define powerLawCurrent_H

#include "waveTheory.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

/*---------------------------------------------------------------------------*\
                       Class powerLawCurrent Declaration
\*---------------------------------------------------------------------------*/

class powerLawCurrent
:
    public waveTheory
{

protected:

    // Protected data
        vector Uref_;        // velocity at (and above) the still water level
        scalar depth_;       // water depth [m]
        scalar exponent_;    // profile exponent (default 1/7)
        scalar Tsoft_;
        scalar localSeaLevel_;

    // Protected member functions

        scalar factor(const scalar&) const;

public:

    //- Runtime type information
    TypeName("powerLawCurrent");

    // Constructors

        //- from components
        powerLawCurrent
        (
            const word&,
            const fvMesh& mesh_
        );


    // Destructor

        ~powerLawCurrent()
        {}


    // Member Functions
        virtual Switch conflictTSoftInitialise() const
        {
            if (Foam::mag(Uref_) > SMALL && Tsoft_ > 0)
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

        //- Depth-AVERAGED speed: Uref/(exponent+1); used by GABC-type callers
        //  that want one representative advection speed.
        inline virtual vector currentSpeed()
        {
            return Uref_/(exponent_ + 1.0);
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
