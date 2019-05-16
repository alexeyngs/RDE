
#ifndef rhoRDEThermo_H
#define rhoRDEThermo_H
#include "moleFractions.H"
#include "rhoReactionThermo.H"
namespace Foam
{
class rhoRDEThermo : public rhoReactionThermo
{
public:
    //- Runtime type information
    TypeName("rhoRDEThermo");

    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        rhoRDEThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );
    //- Construct from mesh and phase name
    rhoRDEThermo(const fvMesh&,const word& phaseName);
    //- Standard selection based on fvMesh
    static autoPtr<rhoRDEThermo> New(const fvMesh&,const word& phaseName=word::null);
    //- Destructor
    virtual ~rhoRDEThermo();
    //---------------------------------------------
    //---------------------------------------------
    virtual volScalarField & MolWeight() = 0;
    virtual volScalarField & FMolWeight() = 0;
    virtual volScalarField & Induction() = 0;
    virtual volScalarField & FInduction() = 0;
    virtual void chemistry(scalar CameraLength) = 0;
    virtual void CorrectMinMax() = 0;
    virtual scalar Dcj() = 0;
    virtual void CorrectErrors(bool UseChemistry, volScalarField & rho, surfaceScalarField & Rhof, surfaceScalarField & Pf) = 0;
    virtual scalar GetmolWeight(Foam::label celli, Foam::scalar Ethermodynamical) = 0;
    virtual Foam::scalar Ethermodynamical(Foam::label celli) = 0;
    virtual scalar Ethermodynamical(scalar P, scalar rho, scalar molWeight) = 0;
    virtual Foam::scalar patchMolWeight(const label patchi, const label facei) = 0;
};

} // End namespace Foam

#endif
