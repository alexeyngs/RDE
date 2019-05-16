#ifndef RDEMixtureEquilibrium_H
#define RDEMixtureEquilibrium_H

namespace Foam{
template<class ThermoType>
class RDEMixtureEquilibrium : public RDEMixtureBasic<ThermoType>
{
    //- Disallow default bitwise copy construct
    RDEMixtureEquilibrium(const RDEMixtureEquilibrium &) = delete;
    //- Disallow default bitwise assignment
    void operator = (const RDEMixtureEquilibrium &) = delete;
public:
    typedef ThermoType thermoType;
    void SetMixtureEquilibrium(const bool zero, const Foam::scalar Rho, const Foam::scalar E, const Foam::scalar Induction, const Foam::scalar alpha, const ThermoType & mixture)const;
    scalar MoleSpecieEquilibrium(const word name, const scalar Rho, const scalar E) const; // Равновесная термодинамика от (Rho, E)
    //----------------------------------------------------------------
    //- Construct from dictionary, mesh and phase name
    RDEMixtureEquilibrium(const dictionary&, const fvMesh&, const word&);
    //- Destructor
    virtual ~RDEMixtureEquilibrium(){}
    //- Return the instantiated type name
    static word typeName()
    {
        return "RDEMixtureEquilibrium<" + ThermoType::typeName() + '>';
    }
    const ThermoType & cellMixture(const label celli) const;
    const ThermoType & patchFaceMixture(const label patchi,const label facei) const;
};
} // End namespace Foam

template<class ThermoType>
Foam::RDEMixtureEquilibrium<ThermoType>::RDEMixtureEquilibrium(const dictionary& thermoDict, const fvMesh& mesh, const word& phaseName) : 
Foam::RDEMixtureBasic<ThermoType>(thermoDict, mesh, phaseName)
{
}

extern const Foam::volScalarField * FieldRho;
extern const Foam::volScalarField * FieldRhoE;


template<class ThermoType>
const ThermoType & Foam::RDEMixtureEquilibrium<ThermoType>::cellMixture(const label celli) const
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
    if(FieldRho && FieldRhoE)
    {
        const Foam::scalar Rho = (*FieldRho)[celli];
        const Foam::scalar E = (*FieldRhoE)[celli] / (*FieldRho)[celli];
        const scalar Sigma = mixture.W()*mixture.moles();           // Суммарная масса всей смеси
        const Foam::scalar MolWeightStar = GetMolWeightInduction(this->species(), this->speciesData(), 0.0, this->powerN, Sigma);
        Foam::scalar alpha = (MolWeight - MolWeightStar)/(this->molWeightCJ - MolWeightStar);
        if(alpha < 0.0) alpha = 0.0;
        if(alpha > 1.0) alpha = 1.0;
        SetMixtureEquilibrium(false, Rho, E, Induction, alpha, mixture);
    }else // Когда ещё не поймали указатели(в начале программы), то считаем, что химии нет
    {
        SetMixtureEquilibrium(true, 0.0, 0.0, Induction, 0.0, mixture);
    }
    return mixture;
}

template<class ThermoType>
const ThermoType & Foam::RDEMixtureEquilibrium<ThermoType>::patchFaceMixture(const label patchi, const label facei) const
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
    if(FieldRho && FieldRhoE)
    {
        const Foam::scalar Rho = FieldRho->boundaryField()[patchi][facei];
        const Foam::scalar E = FieldRhoE->boundaryField()[patchi][facei] / FieldRho->boundaryField()[patchi][facei];
        const scalar Sigma = mixture.W()*mixture.moles();           // Суммарная масса всей смеси
        const Foam::scalar MolWeightStar = GetMolWeightInduction(this->species(), this->speciesData(), 0.0, this->powerN, Sigma);
        Foam::scalar alpha = (MolWeight - MolWeightStar)/(this->molWeightCJ - MolWeightStar);
        if(alpha < 0.0) alpha = 0.0;
        if(alpha > 1.0) alpha = 1.0;
        SetMixtureEquilibrium(false, Rho, E, Induction, alpha, mixture);
    }else // Когда ещё не поймали указатели(в начале программы), то считаем, что химии нет
    {
        SetMixtureEquilibrium(true, 0.0, 0.0, Induction, 0.0, mixture);
    }
    return mixture;
}


template<class ThermoType>
void Foam::RDEMixtureEquilibrium<ThermoType>::SetMixtureEquilibrium(const bool zero, const Foam::scalar Rho, const Foam::scalar E, const Foam::scalar Induction, const Foam::scalar alpha, const ThermoType & mixture)const
{
    // zero - начальное состояние, когда alpha=0, и не ничто не зависит от Rho и E
    // Суммарная масса смеси
    const Foam::scalar Sigma = mixture.W()*mixture.moles();
    if(Induction > SMALL)
    {
        mixture.molWeight_ = GetMolWeightInduction(this->species(), this->speciesData(), Induction, this->powerN, Sigma);
        mixture.Hchemical_ = GetHchemical(this->species(), this->speciesData(), Induction, this->powerN, Sigma);
    }else
    {
        const Foam::scalar MolWeightStar = GetMolWeightInduction(this->species(), this->speciesData(), 0.0, this->powerN, Sigma);
        const Foam::scalar HchemicalStar = GetHchemical(this->species(), this->speciesData(), 0.0, this->powerN, Sigma);
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        mixture.molWeight_ = 0.0;       // молярная масса смеси ячейки
        mixture.Edissociation_ = 0.0;   // энергия диссоцияции ячейки
        Foam::scalar SummaMoles = 0.0;  // суммарное количество молей ячейки
        forAll(this->species(), i)
        {
            Foam::scalar mole = zero ? this->speciesData()[i].moles() : 
                (1.0-alpha)*this->speciesData()[i].moles() + alpha*MoleSpecieEquilibrium(this->species()[i], Rho, E);
            SummaMoles += mole;  // Сумма всех молей
            mixture.molWeight_ += this->speciesData()[i].molWeight_ * mole;
            mixture.Edissociation_ += this->speciesData()[i].Edissociation_ * mole;
        }
        mixture.molWeight_ /= SummaMoles;
        mixture.Edissociation_ /= SummaMoles;
        //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        Foam::scalar Hdelta = HchemicalStar - mixture.Edissociation_*(1.0/MolWeightStar - 1.0/this->MolWeightMin);
        Foam::scalar Psi = (mixture.molWeight_ < MolWeightStar) ? 1.0 : (mixture.molWeight_ > this->molWeightCJ) ? 0.0 : (1.0/this->molWeightCJ - 1.0/mixture.molWeight_) / (1.0/this->molWeightCJ - 1.0/MolWeightStar);
        mixture.Hchemical_ = mixture.Edissociation_*(1.0/mixture.molWeight_ - 1.0/this->MolWeightMin) + Hdelta*Psi;
    }
}

template<typename ThermoType>
Foam::scalar Foam::RDEMixtureEquilibrium<ThermoType>::MoleSpecieEquilibrium(const Foam::word name, const Foam::scalar Rho, const Foam::scalar E) const
{
    return 0.1;
}
#endif