#ifndef RDEMixtureVasilev_H
#define RDEMixtureVasilev_H
namespace Foam{
template<class ThermoType>
class RDEMixtureVasilev : public RDEMixtureBasic<ThermoType>
{
    //- Disallow default bitwise copy construct
    RDEMixtureVasilev(const RDEMixtureVasilev &) = delete;
    //- Disallow default bitwise assignment
    void operator = (const RDEMixtureVasilev &) = delete;
public:
    typedef ThermoType thermoType;
    void SetMixtureVasilev(const Foam::scalar Induction, Foam::scalar & MolWeight, const ThermoType & mixture)const;
    //----------------------------------------------------------------
    //- Construct from dictionary, mesh and phase name
    RDEMixtureVasilev(const dictionary&, const fvMesh&, const word&);
    //- Destructor
    virtual ~RDEMixtureVasilev(){}
    //- Return the instantiated type name
    static word typeName()
    {
        return "RDEMixtureVasilev<" + ThermoType::typeName() + '>';
    }
    const ThermoType & cellMixture(const label celli) const;
    const ThermoType & patchFaceMixture(const label patchi,const label facei) const;
};
} // End namespace Foam

template<class ThermoType>
Foam::RDEMixtureVasilev<ThermoType>::RDEMixtureVasilev(const dictionary& thermoDict, const fvMesh& mesh, const word& phaseName) : 
Foam::RDEMixtureBasic<ThermoType>(thermoDict, mesh, phaseName)
{
}


template<class ThermoType>
const ThermoType & Foam::RDEMixtureVasilev<ThermoType>::cellMixture(const label celli) const
{
    const Foam::IMixture * Mixture = dynamic_cast<const Foam::IMixture*>(this);
    // mixture - одна универсальная переменная
    const ThermoType & mixture = Foam::multiComponentMixture<ThermoType>::cellMixture(celli);
    // MolWeight - молекулярная масса текущей ячейки
    Foam::scalar & MolWeight = ((Foam::scalarField &)this->MolWeight_.primitiveField())[celli];
    // Induction - индукция масса текущей ячейки
    const Foam::scalar Induction = this->Induction_[celli];
    mixture.Mixture = const_cast<Foam::IMixture*>(Mixture);
    mixture.Einduction_ = this->Einduction;
    mixture.Kinduction_ = this->Kinduction;
    // Получение необходимого - Edissociation - molWeight - Hchemical
    mixture.Edissociation_ = this->Edissociation;
    SetMixtureVasilev(Induction, MolWeight, mixture);
    return mixture;
}

template<class ThermoType>
const ThermoType & Foam::RDEMixtureVasilev<ThermoType>::patchFaceMixture(const label patchi, const label facei) const
{
    const Foam::IMixture * Mixture = dynamic_cast<const Foam::IMixture*>(this);
    // mixture - одна универсальная переменная
    const ThermoType & mixture = Foam::multiComponentMixture<ThermoType>::patchFaceMixture(patchi, facei);
    // MolWeight - молекулярная масса текущей ячейки
    Foam::scalar & MolWeight = ((Foam::volScalarField::Boundary  &)this->MolWeight_.boundaryField())[patchi][facei];
    // Induction - индукция масса текущей ячейки
    const Foam::scalar Induction = this->Induction_.boundaryField()[patchi][facei];
    mixture.Mixture = const_cast<Foam::IMixture*>(Mixture);
    mixture.Einduction_ = this->Einduction;
    mixture.Kinduction_ = this->Kinduction;
    // Получение необходимого - Edissociation - molWeight - Hchemical
    mixture.Edissociation_ = this->Edissociation;
    SetMixtureVasilev(Induction, MolWeight, mixture);
    return mixture;
}

template<class ThermoType>
void Foam::RDEMixtureVasilev<ThermoType>::SetMixtureVasilev(const Foam::scalar Induction, Foam::scalar & MolWeight, const ThermoType & mixture)const
{
    // Суммарная масса смеси
    const Foam::scalar Sigma = mixture.W()*mixture.moles();
    if(Induction > SMALL)
    {
        MolWeight = GetMolWeightInduction(this->species(), this->speciesData(), Induction, this->powerN, Sigma);
        mixture.Hchemical_ = GetHchemical(this->species(), this->speciesData(), Induction, this->powerN, Sigma);
    }else
    {
        register Foam::scalar MolWeightStar = GetMolWeightInduction(this->species(), this->speciesData(), 0.0, this->powerN, Sigma);
        register Foam::scalar HchemicalStar = GetHchemical(this->species(), this->speciesData(), 0.0, this->powerN, Sigma);
        register Foam::scalar Hdelta = HchemicalStar - mixture.Edissociation_*(1.0/MolWeightStar - 1.0/this->MolWeightMin);
        register Foam::scalar Psi = (MolWeight < MolWeightStar) ? 1.0 : (MolWeight > this->molWeightCJ) ? 0.0 : (1.0/this->molWeightCJ - 1.0/MolWeight) / (1.0/this->molWeightCJ - 1.0/MolWeightStar);
        mixture.Hchemical_ = mixture.Edissociation_*(1.0/MolWeight - 1.0/this->MolWeightMin) + Hdelta*Psi;
    }
    mixture.molWeight_ = MolWeight;
}
#endif