#ifndef RDEMixtureNikolaev_H
#define RDEMixtureNikolaev_H
namespace Foam{
template<class ThermoType>
class RDEMixtureNikolaev : public RDEMixtureBasic<ThermoType>
{
    //- Disallow default bitwise copy construct
    RDEMixtureNikolaev(const RDEMixtureNikolaev &) = delete;
    //- Disallow default bitwise assignment
    void operator = (const RDEMixtureNikolaev &) = delete;
public:
    typedef ThermoType thermoType;
    //----------------------------------------------------------------
    //- Construct from dictionary, mesh and phase name
    RDEMixtureNikolaev(const dictionary&, const fvMesh&, const word&);
    //- Destructor
    virtual ~RDEMixtureNikolaev(){}
    //- Return the instantiated type name
    static word typeName()
    {
        return "RDEMixtureNikolaev<" + ThermoType::typeName() + '>';
    }
    const ThermoType & cellMixture(const label celli) const;
    const ThermoType & patchFaceMixture(const label patchi,const label facei) const;
};
} // End namespace Foam

template<class ThermoType>
Foam::RDEMixtureNikolaev<ThermoType>::RDEMixtureNikolaev(const dictionary& thermoDict, const fvMesh& mesh, const word& phaseName) : 
Foam::RDEMixtureBasic<ThermoType>(thermoDict, mesh, phaseName)
{
}

template<class ThermoType>
const ThermoType & Foam::RDEMixtureNikolaev<ThermoType>::cellMixture(const label celli) const
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
    mixture.molWeight_ = MolWeight;
    mixture.Hchemical_ = mixture.Edissociation_*(1.0/mixture.molWeight_ - 1.0/this->MolWeightMin);
    return mixture;
}

template<class ThermoType>
const ThermoType & Foam::RDEMixtureNikolaev<ThermoType>::patchFaceMixture(const label patchi, const label facei) const
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
    mixture.molWeight_ = MolWeight;
    mixture.Hchemical_ = mixture.Edissociation_*(1.0/mixture.molWeight_ - 1.0/this->MolWeightMin);
    return mixture;
}
#endif