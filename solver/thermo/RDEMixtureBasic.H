#ifndef RDEMixtureBasic_H
#define RDEMixtureBasic_H

#include "fvCFD.H"
#include "fvMesh.H"
#include "reactingMixture.H"
#include "rhoReactionThermo.H"
#include "../interface.H"
//----------------------------------------------------------------
const scalar Kplus = 6.0e8;                     // K+
const scalar T0 = 300.0;                        // T0
//----------------------------------------------------------------
namespace Foam{
template<class ThermoType>
class RDEMixtureBasic : public IMixture, public speciesTable, public multiComponentMixture<ThermoType>
{
    //- Disallow default bitwise copy construct
    RDEMixtureBasic(const RDEMixtureBasic &) = delete;
    //- Disallow default bitwise assignment
    void operator = (const RDEMixtureBasic &) = delete;
public:
    //- The type of thermo package this mixture is instantiated for
    typedef ThermoType thermoType;
    //- Table of species composition
    speciesCompositionTable speciesComposition_;
    volScalarField MolWeight_;
    volScalarField FMolWeight_;
    volScalarField Induction_;
    volScalarField FInduction_;
    scalar powerN;
    scalar theta;
    scalar beta;
    scalar Einduction;    // Энергия активации [Дж/моль]
    scalar Kinduction;     // Константа активации в моль*с/м3
    scalar molWeightCJ;
    scalar Kminus;
    scalar StartInduction;
    scalar Dcj;
    //----------------------------------------------------------------
    scalar MolWeightA;          // Атомарная масса
    scalar MolWeightMin;        // Масса предельной диссоциации
    scalar MolWeightMax;
    scalar Edissociation;
    //----------------------------------------------------------------
    //- Construct from dictionary, mesh and phase name
    RDEMixtureBasic(const dictionary&, const fvMesh&, const word&);
    //- Destructor
    virtual ~RDEMixtureBasic(){}
    //- Return the instantiated type name
    static word typeName()
    {
        return "RDEMixtureBasic<" + ThermoType::typeName() + '>';
    }
    //- Read dictionary
    void read(const dictionary&);
    //- Table of species composition
    const speciesCompositionTable & specieComposition() const
    {
        return speciesComposition_;
    }
    Foam::scalar W(const label speciei) const
    {
        const ThermoType & type = multiComponentMixture<ThermoType>::getLocalThermo(speciei);
        return type.W()/type.Y();
    }
    Foam::scalar GetMolWeight(Foam::label celli, Foam::scalar Hchemical)const;
    void CorrectMinMax();
    std::pair<scalar, scalar> SetWeightAndT(const ThermoType & mixture, const scalar h, const scalar rho, const scalar Pav, const scalar molWeightOld) const;
    void Chemistry(const label celli, const scalar rho, const scalar T);
    //-------------------------------------------------
    scalar GetA(const scalar MolWeight, const scalar T) const override final;
    scalar Cv(const scalar P, const scalar T, const scalar MolWeight)const override final;
    scalar GetTmax(const scalar MolWeight, const scalar Uthermodinamical) const override final;
};
} // End namespace Foam
//---------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------
#include "VasilevChemistry.H"


template<class ThermoType>
Foam::RDEMixtureBasic<ThermoType>::RDEMixtureBasic(const dictionary& thermoDict, const fvMesh& mesh, const word& phaseName) : 
    speciesTable(),
    multiComponentMixture<ThermoType>
    (
        thermoDict,
        mesh,
        phaseName
    ),
    speciesComposition_(),
    MolWeight_
    (
        IOobject
        (
            "MolWeight",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    FMolWeight_
    (
        IOobject
        (
            "FMolWeight",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("dimVelocity*dimMoles", dimensionSet(1, -3, -1, 0, 1, 0, 0), 0.0)
    ),
    Induction_
    (
        IOobject
        (
            "Induction",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    FInduction_
    (
        IOobject
        (
            "FInduction",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("dimVelocity", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
    )
{
    const Foam::Time & runTime = mesh.time();
    const ThermoType & mixture = multiComponentMixture<ThermoType>::cellMixture(0);
    //-----------------------------------------------------------------
    MolWeightA = mixture.W()/mixture.atoms();                   // Общее количество атомов с учетом молей
    MolWeightMin = mixture.W()/mixture.decomposition();         // Масса предельной диссоциации
    const scalar fuel = mixture.needO()*mixture.moles();        // Количество атомов кислорода для полного сгорания
    const scalar Sigma = mixture.W()*mixture.moles();           // Суммарная масса всей смеси
    const scalar atoms = mixture.atoms()*mixture.moles();       // Общее количество атомов с учетом молей
    // Получение параметров, которые зависят от состава смеси
    powerN = readScalar(thermoDict.subDict("thermo").lookup("powerN"));
    scalar MolWeight0 = GetMolWeightInduction<ThermoType>(this->species(), this->speciesData(), 1.0, powerN, Sigma);
    scalar MolWeightStar = GetMolWeightInduction<ThermoType>(this->species(), this->speciesData(), 0.0, powerN, Sigma);
    MolWeightMax = GetMolWeightMax<ThermoType>(this->species(), this->speciesData(), fuel, Sigma);
    Edissociation = GetEdissociation<ThermoType>(this->species(), this->speciesData(), fuel, Sigma, MolWeightMin);
    //-----------------------------------------------------------------
    Info << "Edissociation = " << Edissociation << " J/mol (" << Edissociation/4186.8 << " kcal/mol)" << ", Needs atoms(Oxigen) to total dissociate = " << mixture.needO()*mixture.moles() << endl;
    Info << "Sigma = " << Sigma << ", MolWeight0 = " << MolWeight0 << ", MolWeightStar = " << MolWeightStar << ", MolWeightA = " << MolWeightA << ", MolWeightMin = " << MolWeightMin << ", MolWeightMax = " << MolWeightMax << endl;
    forAll(multiComponentMixture<ThermoType>::speciesData(), i)
    {
        multiComponentMixture<ThermoType>::setActive(i);
    }

    // Заполнение поля молярной массы
    if(runTime.value() < Foam::SMALL) 
    {
        // Заполнение поля в центре области
        forAll(basicMultiComponentMixture::Y(0) , i)
        {
            const ThermoType & mixture = multiComponentMixture<ThermoType>::cellMixture(i);
            scalar Induction = Induction_[i];
            MolWeight_[i] = GetMolWeightInduction<ThermoType>(this->species(), this->speciesData(), Induction, powerN, mixture.W()*mixture.moles());
        }
        // Заполнение поля на границах
        forAll(MolWeight_.boundaryField(), patchi) forAll(MolWeight_.boundaryField()[patchi], facei)
        {
            MolWeight_.boundaryFieldRef()[patchi][facei] = MolWeight0;
        }
    }

    theta = readScalar(thermoDict.subDict("thermo").lookup("theta"));
    beta = readScalar(thermoDict.subDict("thermo").lookup("beta"));
    molWeightCJ = readScalar(thermoDict.subDict("thermo").lookup("molWeightCJ"));
    Kminus = readScalar(thermoDict.subDict("thermo").lookup("Kminus"));
    StartInduction = readScalar(thermoDict.subDict("thermo").lookup("StartInduction"));
    Dcj = readScalar(thermoDict.subDict("thermo").lookup("Dcj"));
    Edissociation *= 1000.0;
    //----------------------------------------------------------------
    // Индукция
    scalar SumMoles = 0.0;
    Kinduction = 1.0;
    forAll(this->speciesData(), i) if(this->speciesData()[i].Kinduction_ > 0.0) SumMoles += this->speciesData()[i].moles();
    forAll(this->speciesData(), i) if(this->speciesData()[i].Kinduction_ > 0.0)
    {
        scalar beta = this->speciesData()[i].moles() / SumMoles;
        Kinduction *= Foam::pow(this->speciesData()[i].Kinduction_, beta);
    }
    SumMoles = 0.0;
    Einduction = 0.0;
    forAll(this->speciesData(), i) if(this->speciesData()[i].Einduction_ > 0.0) SumMoles += this->speciesData()[i].moles();
    forAll(this->speciesData(), i) if(this->speciesData()[i].Einduction_ > 0.0)
    {
        scalar beta = this->speciesData()[i].moles() / SumMoles;
        Einduction += this->speciesData()[i].Einduction_ * beta;
    }
    Info << "Einduction = " << Einduction << " J/mol, Kinduction = " << Kinduction << " mol*c/m3" << endl;
    Info << "(none)Einduction = " << mixture.Einduction() << " J/mol, Kinduction = " << mixture.Kinduction() << " mol*c/m3" << endl;
    Einduction *= 1000.0;
}

template<class ThermoType>
void Foam::RDEMixtureBasic<ThermoType>::read(const dictionary & thermoDict)
{
}

template<class ThermoType>
Foam::scalar Foam::RDEMixtureBasic<ThermoType>::GetMolWeight(Foam::label celli, Foam::scalar Hchemical)const
{
    return 1.0/(1.0/MolWeightMin + Hchemical/Edissociation);
}

template<class ThermoType>
void Foam::RDEMixtureBasic<ThermoType>::CorrectMinMax()
{
    forAll(MolWeight_, celli)
    {
        if(MolWeight_[celli] < MolWeightMin) MolWeight_[celli] = MolWeightMin;
        if(MolWeight_[celli] > MolWeightMax) MolWeight_[celli] = MolWeightMax;
    }
}

template<class ThermoType>
inline std::pair<scalar, scalar> Foam::RDEMixtureBasic<ThermoType>::SetWeightAndT(const ThermoType & mixture, const scalar h, const scalar rho, const scalar Pav, const scalar molWeightOld)const
{
    Foam::scalar Tav = Pav*molWeightOld/rho/Foam::constant::thermodynamic::RR;
    Foam::scalar Uthermodinamical = Pav/rho * this->GetA(molWeightOld, Tav);
    Foam::scalar Hchemistry = h - Uthermodinamical;
    //--------------------------------
    Foam::scalar molWeight = 1.0/(1.0/MolWeightMin + Hchemistry/Edissociation);
    if(molWeight < MolWeightMin) molWeight = MolWeightMin;
    if(molWeight > MolWeightMax) molWeight = MolWeightMax;
    // далее идет GetT ----------------------------------------
    Foam::scalar a,b,U,T;
    a = 50.0;
    b = GetTmax(Uthermodinamical, molWeight);
    while(abs(a-b)/(a+b) > SMALL)
    {
	    T = (a+b)/2.0;
	    U = Foam::constant::thermodynamic::RR*T/molWeight * this->GetA(molWeight, T);
	    if(U > Uthermodinamical) b = T;
        else a = T;
    }
    T = (a+b)/2.0;
    //---------------------------------------------------------
    Info << ", Uthermodinamical = " << Uthermodinamical;
    return std::pair<scalar, scalar> (molWeight, T);
}

template<class ThermoType>
void Foam::RDEMixtureBasic<ThermoType>::Chemistry(const label celli, const scalar rho, const scalar T)
{
    scalar & Induction = Induction_[celli];
    scalar & FInduction = FInduction_[celli];
    scalar & MolWeight = MolWeight_[celli];
    scalar & FMolWeight = FMolWeight_[celli];
    const ThermoType & mixture = Foam::RDEMixtureBasic<ThermoType>::cellMixture(celli);
    const scalar fuel = mixture.needO()*mixture.moles();                         // Количество атомов кислорода для полного сгорания
    const scalar Sigma = mixture.W()*mixture.moles();                            // Суммарная масса всей смеси
    const scalar atoms = mixture.atoms()*mixture.moles();                        // Общее количество атомов с учетом молей
    //------------------------------------------------------
    // Стадия №1 - задержка воспламенения - уменьшение Induction от 1 до 0
    if(T > StartInduction)
    {
        const label index_O2 = this->species()["O2"];
        if(index_O2 < 0) FatalError << "Не найден индекс кислорода!" << endl << exit(FatalError);
        scalar molesO2 = this->speciesData()[index_O2].moles();
        const label index_H2 = this->species()["H2"];
        scalar molesH2 = this->speciesData()[index_H2].moles();
        register scalar concentrationO2 = rho * molesO2 / Sigma * 1000.0;     // Концентрация кислорода	[кг/м3 * моль / кг] в моль/м3
        register scalar concentrationH2 = rho * molesH2 / Sigma * 1000.0;     // Концентрация водорода	[кг/м3 * моль / кг] в моль/м3
//register scalar tind = Kind / Foam::sqrt(concentrationO2) / Foam::sqrt(concentrationH2) * Foam::exp(Ea / Foam::constant::thermodynamic::RR / T);              // White вычисление время индукции
        register scalar tind = Kinduction / concentrationO2 * Foam::exp(Einduction / Foam::constant::thermodynamic::RR / T);   // Strelow вычисление время индукции
        FInduction = -1.0/tind;
    }else
    {
        FInduction = 0.0;
    }

    // Стадия №2 - химическая реакция
    if(Induction > Foam::SMALL || T < 300.0)
    {
        FMolWeight = 0.0;
    }
    else    // Вычисление приращения молекулярной массы
    {
        register scalar FWeight =
            4.0*Kplus*rho*rho/MolWeight * (1.0 - MolWeight/MolWeightMax)*(1.0 - MolWeight/MolWeightMax) - 
            4.0*Kplus*Kminus*rho*(MolWeight/MolWeightMin - 1.0)*Foam::pow(T/T0, beta/2.0)*Foam::pow(1.0 - Foam::exp(-theta/T), beta)*Foam::exp(-Edissociation/Foam::constant::thermodynamic::RR/T);
		if(FWeight < -20.0) FWeight = -20.0;
		FMolWeight = rho*FWeight;
    }
}

template<class ThermoType>
Foam::scalar Foam::RDEMixtureBasic<ThermoType>::GetA(const Foam::scalar MolWeight, const Foam::scalar T) const
{
    if(T <= 50.0) return 3.0/4.0*(MolWeight/MolWeightA + 1.0); // Если слишком малое T, => огромное theta/T
	return 3.0/4.0*(MolWeight/MolWeightA + 1.0) + 3.0/2.0*(MolWeight/MolWeightA - 1.0)*(theta/T)/(Foam::exp(theta/T) - 1.0);            // Николаев
}

template<class ThermoType>
Foam::scalar Foam::RDEMixtureBasic<ThermoType>::Cv(const Foam::scalar P, const Foam::scalar T, const Foam::scalar MolWeight)const
{
    const scalar R = Foam::constant::thermodynamic::RR / MolWeight;
    return this->GetA(MolWeight, T) * R;
}

template<class ThermoType>
Foam::scalar Foam::RDEMixtureBasic<ThermoType>::GetTmax(const Foam::scalar MolWeight, const Foam::scalar Uthermodinamical) const
{
    return 5000.0;
    Foam::scalar A = 3.0/4.0*(MolWeight/MolWeightA + 1.0);
    return Uthermodinamical / A;
}
#endif
