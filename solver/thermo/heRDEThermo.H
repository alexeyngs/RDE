#ifndef heRDEThermo_H
#define heRDEThermo_H

#include "rhoRDEThermo.H"
#include "RDEThermo.H"
#include "turbulentFluidThermoModel.H"

namespace Foam
{
template<class BasicPsiThermo, class MixtureType>
class heRDEThermo : public heThermo<BasicPsiThermo, MixtureType>
{
    //- Calculate the thermo variables
    void calculate();
    //- Construct as copy (not implemented)
    heRDEThermo(const heRDEThermo<BasicPsiThermo, MixtureType>&);
public:
    //- Runtime type information
    TypeName("heRDEThermo");
    // Constructors
        //- Construct from mesh and phase name
        heRDEThermo(const fvMesh&,const word& phaseName);
    //- Destructor
    virtual ~heRDEThermo();
    // Member functions
    virtual void correct();
    // от класса rhoRDEThermo
    virtual volScalarField & MolWeight() override;
    virtual volScalarField & FMolWeight() override;
    virtual volScalarField & Induction() override;
    virtual volScalarField & FInduction() override;
    virtual void chemistry(scalar CameraLength) override;
    virtual void CorrectMinMax() override;
    virtual Foam::scalar Dcj() override;
    virtual void CorrectErrors(bool UseChemistry, volScalarField & rho, surfaceScalarField & Rhof, surfaceScalarField & Pf) override;
    virtual Foam::scalar GetmolWeight(Foam::label celli, Foam::scalar Ethermodynamical) override;
    virtual Foam::scalar Ethermodynamical(Foam::label celli) override;
    virtual Foam::scalar Ethermodynamical(Foam::scalar P, Foam::scalar rho, Foam::scalar molWeight) override;
    virtual Foam::scalar patchMolWeight(const label patchi, const label facei) override;
};
} // End namespace Foam

template<class BasicPsiThermo, class MixtureType>
void Foam::heRDEThermo<BasicPsiThermo, MixtureType>::calculate()
{
    const scalarField & hCells = this->he();
    const scalarField & pCells = this->p_;
    scalarField & TCells = this->T_.primitiveFieldRef();
    scalarField & psiCells = this->psi_.primitiveFieldRef();
    scalarField & rhoCells = this->rho_.primitiveFieldRef();
    scalarField & muCells = this->mu_.primitiveFieldRef();
    scalarField & alphaCells = this->alpha_.primitiveFieldRef();
    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType & mixture_ = this->cellMixture(celli);
        TCells[celli] = mixture_.THE(hCells[celli], pCells[celli], TCells[celli]);
        psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]);
        rhoCells[celli] = mixture_.rho(pCells[celli], TCells[celli]);
        muCells[celli] = mixture_.mu(pCells[celli], TCells[celli]);
        alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);
    }

    // Граничные условия
    volScalarField::Boundary & pBoundaryField = this->p_.boundaryFieldRef();
    volScalarField::Boundary & TBoundaryField = this->T_.boundaryFieldRef();
    volScalarField::Boundary & psiBoundaryField = this->psi_.boundaryFieldRef();
    volScalarField::Boundary & rhoBoundaryField = this->rho_.boundaryFieldRef();
    volScalarField::Boundary & heBoundaryField = this->he().boundaryFieldRef();
    volScalarField::Boundary & muBoundaryField = this->mu_.boundaryFieldRef();
    volScalarField::Boundary & alphaBoundaryField = this->alpha_.boundaryFieldRef();
    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField & pp = pBoundaryField[patchi];
        fvPatchScalarField & pT = TBoundaryField[patchi];
        fvPatchScalarField & ppsi = psiBoundaryField[patchi];
        fvPatchScalarField & prho = rhoBoundaryField[patchi];
        fvPatchScalarField & phe = heBoundaryField[patchi];
        fvPatchScalarField & pmu = muBoundaryField[patchi];
        fvPatchScalarField & palpha = alphaBoundaryField[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType & mixture_ = this->patchFaceMixture(patchi, facei);
                phe[facei] = mixture_.HE(pp[facei], pT[facei]);
                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                prho[facei] = mixture_.rho(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType & mixture_ = this->patchFaceMixture(patchi, facei);
                pT[facei] = mixture_.THE(phe[facei], pp[facei], pT[facei]);
                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                prho[facei] = mixture_.rho(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::heRDEThermo<BasicPsiThermo, MixtureType>::heRDEThermo(const fvMesh& mesh, const word& phaseName):
heThermo<BasicPsiThermo, MixtureType>(mesh, phaseName)
{
    calculate();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::heRDEThermo<BasicPsiThermo, MixtureType>::~heRDEThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::heRDEThermo<BasicPsiThermo, MixtureType>::correct()
{
    if (debug)
    {
        InfoInFunction << endl;
    }
    calculate();
    if (debug)
    {
        Info<< "    Finished" << endl;
    }
    
}

template<class BasicPsiThermo, class MixtureType>
Foam::volScalarField & Foam::heRDEThermo<BasicPsiThermo, MixtureType>::MolWeight()
{
    return MixtureType::MolWeight_;
}

template<class BasicPsiThermo, class MixtureType>
Foam::volScalarField & Foam::heRDEThermo<BasicPsiThermo, MixtureType>::FMolWeight()
{
    return MixtureType::FMolWeight_;
}

template<class BasicPsiThermo, class MixtureType>
Foam::volScalarField & Foam::heRDEThermo<BasicPsiThermo, MixtureType>::Induction()
{
    return MixtureType::Induction_;
}

template<class BasicPsiThermo, class MixtureType>
Foam::volScalarField & Foam::heRDEThermo<BasicPsiThermo, MixtureType>::FInduction()
{
    return MixtureType::FInduction_;
}

template<class BasicPsiThermo, class MixtureType>
void Foam::heRDEThermo<BasicPsiThermo, MixtureType>::chemistry(scalar CameraLength)
{ 
    const scalarField & TCells = this->T_.primitiveField();
    const scalarField & rhoCells = this->rho_.primitiveField();
    const fvMesh & mesh = this->T_.mesh();
    forAll(TCells, celli)
    {
        const point & Position = mesh.C()[celli];
        if(Position.z() < CameraLength)
            MixtureType::Chemistry(celli, rhoCells[celli], TCells[celli]);
    }
}

template<class BasicPsiThermo, class MixtureType>
void Foam::heRDEThermo<BasicPsiThermo, MixtureType>::CorrectMinMax()
{
    MixtureType::CorrectMinMax();
}

template<class BasicPsiThermo, class MixtureType>
Foam::scalar Foam::heRDEThermo<BasicPsiThermo, MixtureType>::Dcj()
{
    return MixtureType::Dcj;
}

template<class BasicPsiThermo, class MixtureType>
void Foam::heRDEThermo<BasicPsiThermo, MixtureType>::CorrectErrors(bool UseChemistry, volScalarField & rho, surfaceScalarField & Rhof, surfaceScalarField & Pf)
{
    const scalarField & hCells = this->he();
    scalarField & rhoCells = this->rho_.primitiveFieldRef();
    scalarField & TCells = this->T_.primitiveFieldRef();
    scalarField & MolWeight = MixtureType::MolWeight_;
    scalarField & pCells = this->p_.primitiveFieldRef();
    scalarField & psiCells = this->psi_.primitiveFieldRef();
    tmp<volScalarField> Pav = fvc::average(Pf);
    tmp<volScalarField> Rhoav = fvc::average(Rhof);
    // есть варианты fvc::surfaceSum(Pf), fvc::reconstruct(Pf)
    forAll(TCells, celli)
    {
        if(TCells[celli] < 50.0 || rho[celli] < 0.0)
        {
            Info << "T < 50K!(celli = " << celli << ")Pold = " << pCells[celli]/101325.0 << ", Told = " << TCells[celli] <<
            ", MolWeightold = " << MolWeight[celli];
            //-----------------------------------------------------------
            pCells[celli] = Pav.ref()[celli];
            if(rho[celli] < 0.0)
            {
                rho[celli] = Rhoav.ref()[celli];
                Info << "(ERROR! RHO < 0) ";
                if(rho[celli] < 0.0)
                {
                    Info << "ERROR! RHO < 0 AGAIN!" << endl;
                    rho[celli] = rho[celli+1];
                }
            }
            // установим молярную массу и температуру
            const typename MixtureType::thermoType & mixture = this->cellMixture(celli);
            std::pair<scalar, scalar> WeightAndT = MixtureType::SetWeightAndT(mixture, hCells[celli],  rho[celli], pCells[celli], MolWeight[celli]);
            MolWeight[celli] = WeightAndT.first;
            TCells[celli] = WeightAndT.second;
            //-------------------------------------------------------------
            psiCells[celli] = mixture.psi(pCells[celli], TCells[celli]);
            rhoCells[celli] = mixture.rho(pCells[celli], TCells[celli]);
            Info << ", P = " << pCells[celli]/101325.0 << ", rho = " << rho[celli] << 
            ", MolWeight = " << MolWeight[celli] << ", T = " << TCells[celli] << endl;
        }
    }
    Pav.clear();
    Rhoav.clear();

    // Граничные условия
    volScalarField::Boundary & pBf = this->p_.boundaryFieldRef();
    volScalarField::Boundary & TBf = this->T_.boundaryFieldRef();
    volScalarField::Boundary & MolWeightBf = MixtureType::MolWeight_.boundaryFieldRef();
    volScalarField::Boundary & psiBf = this->psi_.boundaryFieldRef();
    volScalarField::Boundary & rhoBf = this->rho_.boundaryFieldRef();
    volScalarField::Boundary & heBf = this->he().boundaryFieldRef();
    volScalarField::Boundary & muBf = this->mu_.boundaryFieldRef();
    volScalarField::Boundary & alphaBf = this->alpha_.boundaryFieldRef();
    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField & pp = pBf[patchi];
        fvPatchScalarField & pT = TBf[patchi];
        fvPatchScalarField & pMolWeight = MolWeightBf[patchi];
        fvPatchScalarField & ppsi = psiBf[patchi];
        fvPatchScalarField & prho = rhoBf[patchi];
        fvPatchScalarField & phe = heBf[patchi];
        fvPatchScalarField & pmu = muBf[patchi];
        fvPatchScalarField & palpha = alphaBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei) if(pT[facei] < 50.0 || prho[facei] < 0.0)
            {
                Info << "T < 50K!; BOUNDS: patchi = " << patchi << ", facei = " << facei <<
                " Pold = " << pp[facei]/101325.0 << ", Told = " << pT[facei] << ", MolWeightold = " << pMolWeight[facei];
                //-----------------------------------------------------------
                // установим молярную массу и температуру
                const typename MixtureType::thermoType & mixture = this->patchFaceMixture(patchi, facei);
                std::pair<scalar, scalar> WeightAndT = MixtureType::SetWeightAndT(mixture, phe[facei], prho[facei], pp[facei], pMolWeight[facei]);
                pMolWeight[facei] = WeightAndT.first;
                pT[facei] = WeightAndT.second;
                //-------------------------------------------------------------
                ppsi[facei] = mixture.psi(pp[facei], pT[facei]);
                prho[facei] = mixture.rho(pp[facei], pT[facei]);
                Info << ", P = " << pp[facei]/101325.0 << ", rho = " << prho[facei] << 
                ", MolWeight = " << pMolWeight[facei] << ", T = " << pT[facei] << endl;
            }
        }
        else
        {
            forAll(pT, facei) if(pT[facei] < 50.0 || prho[facei] < 0.0)
            {
                Info << "T < 50K!; BOUNDS: patchi = " << patchi << ", facei = " << facei <<
                " Pold = " << pp[facei]/101325.0 << ", Told = " << pT[facei] << ", MolWeightold = " << pMolWeight[facei];
                //-----------------------------------------------------------
                // установим молярную массу и температуру
                const typename MixtureType::thermoType & mixture = this->patchFaceMixture(patchi, facei);
                std::pair<scalar, scalar> WeightAndT = MixtureType::SetWeightAndT(mixture, phe[facei], prho[facei], pp[facei], pMolWeight[facei]);
                pMolWeight[facei] = WeightAndT.first;
                pT[facei] = WeightAndT.second;
                //pT[facei] = MixtureType::SetWeightAndT();
                //-------------------------------------------------------------
                ppsi[facei] = mixture.psi(pp[facei], pT[facei]);
                prho[facei] = mixture.rho(pp[facei], pT[facei]);
                Info << ", P = " << pp[facei]/101325.0 << ", rho = " << prho[facei] << 
                ", MolWeight = " << pMolWeight[facei] << ", T = " << pT[facei] << endl;
            }
        }
    }
    //--------------------------------------------------

    forAll(TCells, celli) if(TCells[celli] < 50.0 || pCells[celli] < 0.0 || rho[celli] < 0.0)
    {
        scalar T = TCells[celli];
        scalar P = pCells[celli];
        scalar rhorho = rho[celli];
        Foam::Info << "ERROR! AGAIN! celli = " << celli << " T = " << T << " , P = " << P/101325.0 << " , rho = " << rhorho << Foam::endl;
    }
}

template<class BasicPsiThermo, class MixtureType>
Foam::scalar Foam::heRDEThermo<BasicPsiThermo, MixtureType>::GetmolWeight(Foam::label celli, Foam::scalar Ethermodynamical)
{
    const scalarField & hCells = this->he();
    return MixtureType::GetMolWeight(celli, hCells[celli] - Ethermodynamical);
}

template<class BasicPsiThermo, class MixtureType>
Foam::scalar Foam::heRDEThermo<BasicPsiThermo, MixtureType>::Ethermodynamical(Foam::scalar P, Foam::scalar rho, Foam::scalar MolWeight)
{
    Foam::scalar T = P*MolWeight/rho/Foam::constant::thermodynamic::RR;
    const typename MixtureType::thermoType & mixture = this->cellMixture(0);
    mixture.molWeight_ = MolWeight;
    return mixture.Hs(P,T);
}

template<class BasicPsiThermo, class MixtureType>
Foam::scalar Foam::heRDEThermo<BasicPsiThermo, MixtureType>::Ethermodynamical(Foam::label celli)
{
    const scalarField & hCells = this->he();
    const typename MixtureType::thermoType & mixture = this->cellMixture(celli);
    return hCells[celli] - mixture.Hc();
}

template<class BasicPsiThermo, class MixtureType>
Foam::scalar Foam::heRDEThermo<BasicPsiThermo, MixtureType>::patchMolWeight(const Foam::label patchi, const Foam::label facei)
{
    const typename MixtureType::thermoType & mixture = this->patchFaceMixture(patchi, facei);
    return mixture.W();
}

#endif
