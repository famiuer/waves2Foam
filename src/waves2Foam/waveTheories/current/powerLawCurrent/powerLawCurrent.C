/*---------------------------------------------------------------------------*\
    powerLawCurrent — see powerLawCurrent.H
\*---------------------------------------------------------------------------*/

#include "powerLawCurrent.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(powerLawCurrent, 0);
addToRunTimeSelectionTable(waveTheory, powerLawCurrent, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

powerLawCurrent::powerLawCurrent
(
    const word& subDictName,
    const fvMesh& mesh_
)
:
    waveTheory(subDictName, mesh_),
    Uref_(vector(coeffDict_.lookup("Uref"))),
    depth_(readScalar(coeffDict_.lookup("depth"))),
    exponent_
    (
        coeffDict_.lookupOrDefault<scalar>("exponent", 1.0/7.0)
    ),
    Tsoft_(readScalar(coeffDict_.lookup("Tsoft"))),
    localSeaLevel_
    (
        coeffDict_.lookupOrDefault<scalar>("localSeaLevel", seaLevel_)
    )
{}


void powerLawCurrent::printCoeffs()
{
    Info << "Loading wave theory: " << typeName
         << " (Uref = " << Uref_
         << ", depth = " << depth_
         << ", exponent = " << exponent_ << ")" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar powerLawCurrent::factor(const scalar& time) const
{
    scalar factor(1);
    if (Tsoft_ > 0.0)
    {
        factor = Foam::sin(PI_/2.0/Tsoft_*Foam::min(Tsoft_, time));
    }

    return factor;
}


scalar powerLawCurrent::eta
(
    const point& x,
    const scalar& time
) const
{
    return localSeaLevel_;
}


scalar powerLawCurrent::pExcess
(
    const point& x,
    const scalar& time
) const
{
    return referencePressure(localSeaLevel_);
}


vector powerLawCurrent::U
(
    const point& x,
    const scalar& time
) const
{
    // Normalised elevation above the bed: 0 at z = localSeaLevel - depth,
    // 1 at the still water level; clipped so the profile is Uref above the
    // SWL (air side / wave crests) and 0 below the bed.
    scalar s = (x.z() - (localSeaLevel_ - depth_))/depth_;
    s = Foam::max(scalar(0), Foam::min(scalar(1), s));

    return (Uref_*Foam::pow(s, exponent_)*factor(time));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
