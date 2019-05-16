#include "rhoRDEThermo.H"
#include "fvMesh.H"
#include "specie.H"
namespace Foam
{
    defineTypeNameAndDebug(rhoRDEThermo, 0);
    defineRunTimeSelectionTable(rhoRDEThermo, fvMesh);
}
Foam::rhoRDEThermo::rhoRDEThermo(const fvMesh& mesh,const word& phaseName) : Foam::rhoReactionThermo(mesh, phaseName)
{
}
Foam::autoPtr<Foam::rhoRDEThermo> Foam::rhoRDEThermo::New(const fvMesh& mesh,const word& phaseName)
{
    return basicThermo::New<rhoRDEThermo>(mesh, phaseName);
}
Foam::rhoRDEThermo::~rhoRDEThermo()
{
}
