#ifndef RDEThermo_H
#define RDEThermo_H
#include "../interface.H"
namespace Foam
{
// Forward declaration of friend functions and operators
template<class EquationOfState> class RDEThermo;
template<class EquationOfState> inline RDEThermo<EquationOfState> operator+ (const RDEThermo<EquationOfState>&,const RDEThermo<EquationOfState>&);
template<class EquationOfState> inline RDEThermo<EquationOfState> operator*(const scalar, const RDEThermo<EquationOfState>&);
template<class EquationOfState> inline RDEThermo<EquationOfState> operator==(const RDEThermo<EquationOfState>&,const RDEThermo<EquationOfState>&);
template<class EquationOfState> Ostream& operator<<(Ostream&,const RDEThermo<EquationOfState>&);


// class
template<typename EquationOfState>
class RDEThermo : public EquationOfState
{
    // Private
    //- Construct from components
    inline RDEThermo(const EquationOfState & st, const scalar cp, const scalar hf);
public:
    mutable scalar Hchemical_;
    mutable IMixture * Mixture = nullptr;
    //------------------------------------------
    // Constructors
    //- Construct from dictionary
    RDEThermo(const dictionary& dict);
    inline RDEThermo(const EquationOfState & st);
    //- Construct as named copy
    inline RDEThermo(const word&, const RDEThermo&);
    //- Construct and return a clone
    inline autoPtr<RDEThermo> clone() const;
    //- Selector from dictionary
    inline static autoPtr<RDEThermo> New(const dictionary& dict);

    // Member Functions
    //- Return the instantiated type name
    static word typeName()
    {
        return "RDEThermo<" + EquationOfState::typeName() + '>';
    }
//--------------------------------------------
    //- Limit the temperature to be in the range Tlow_ to Thigh_
    inline scalar limit(const scalar T) const;
    //- Heat capacity at constant pressure [J/(kg K)]
    inline scalar Cp(const scalar p, const scalar T) const;
    //- Absolute Enthalpy [J/kg]
    inline scalar Ha(const scalar p, const scalar T) const;
    //- Sensible enthalpy [J/kg]
    inline scalar Hs(const scalar p, const scalar T) const;
    //- Chemical enthalpy [J/kg]
    inline scalar Hc() const;
    //- Энтропия [J/(kg K)]
    inline scalar S(const scalar p, const scalar T) const;
    //- Write to Ostream
    void write(Ostream & os) const;
    // Member operators
    inline void operator+=(const RDEThermo&);
    // Friend operators
    friend RDEThermo operator + <EquationOfState> (const RDEThermo&, const RDEThermo &);
    friend RDEThermo operator * <EquationOfState> (const scalar, const RDEThermo &);
    friend RDEThermo operator == <EquationOfState> (const RDEThermo &, const RDEThermo &);
    // IOstream Operators
    friend Ostream & operator << <EquationOfState> (Ostream &, const RDEThermo &);
};
} // End namespace Foam

//--------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------
// Constructors
template<class EquationOfState>
Foam::RDEThermo<EquationOfState>::RDEThermo(const dictionary & dict) : EquationOfState(dict)
{
}

template<class EquationOfState>
inline Foam::RDEThermo<EquationOfState>::RDEThermo(const EquationOfState & st): EquationOfState(st)
{
}

template<class EquationOfState>
inline Foam::RDEThermo<EquationOfState>::RDEThermo(const word & name, const RDEThermo & ct) : EquationOfState(name, ct)
{
}

template<class EquationOfState>
inline Foam::autoPtr<Foam::RDEThermo<EquationOfState>> Foam::RDEThermo<EquationOfState>::clone() const
{
    return autoPtr<RDEThermo<EquationOfState>>(new RDEThermo<EquationOfState>(*this));
}

template<class EquationOfState>
inline Foam::autoPtr<Foam::RDEThermo<EquationOfState>> Foam::RDEThermo<EquationOfState>::New(const dictionary& dict)
{
    return autoPtr<RDEThermo<EquationOfState>>(new RDEThermo<EquationOfState>(dict));
}

// Member Functions
template<class EquationOfState>
void Foam::RDEThermo<EquationOfState>::write(Ostream& os) const
{
    EquationOfState::write(os);
//    dictionary dict("thermodynamics");
//    dict.add("Edissociation", Edissociation_);
//    os  << indent << dict.dictName() << dict;
}

template<class EquationOfState>
inline Foam::scalar Foam::RDEThermo<EquationOfState>::limit(const scalar T) const
{
    scalar TT = T;
    if(TT < 40.0) TT = 40.0;
    return T;
}

template<class EquationOfState>
inline Foam::scalar Foam::RDEThermo<EquationOfState>::Cp(const scalar p, const scalar T) const
{
    scalar MolWeight = this->W();
    return EquationOfState::CpMCv(p, T) + this->Mixture->Cv(p, T, MolWeight);
}

template<class EquationOfState>
inline Foam::scalar Foam::RDEThermo<EquationOfState>::Ha(const scalar p, const scalar T) const
{
    return Hs(p,T) + Hc();
}

template<class EquationOfState>
inline Foam::scalar Foam::RDEThermo<EquationOfState>::Hs(const scalar p, const scalar T) const
{
    const Foam::scalar MolWeight = this->W();
    return p/EquationOfState::rho(p,T) * this->Mixture->GetA(MolWeight, T);
}

template<class EquationOfState>
inline Foam::scalar Foam::RDEThermo<EquationOfState>::Hc() const
{
    return Hchemical_;
}

template<class EquationOfState>
inline Foam::scalar Foam::RDEThermo<EquationOfState>::S(const scalar p, const scalar T) const
{
    return EquationOfState::S(p, T);
}

//======================================================================================================================================

// Operators  
template<class EquationOfState> 
Foam::Ostream& Foam::operator << (Ostream& os,const RDEThermo<EquationOfState>& ct)
{
    ct.write(os);
    return os;
}

template<class EquationOfState> inline void Foam::RDEThermo<EquationOfState>::operator+=(const RDEThermo<EquationOfState>& ct)
{
    scalar Y1 = this->Y();

    EquationOfState::operator+=(ct);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        scalar Y2 = ct.Y()/this->Y();
    }
}

// Friend Operators
template<class EquationOfState>
inline Foam::RDEThermo<EquationOfState> Foam::operator+(const RDEThermo<EquationOfState> & ct1, const RDEThermo<EquationOfState> & ct2)
{
    EquationOfState eofs(static_cast<const EquationOfState&>(ct1) + static_cast<const EquationOfState&>(ct2));
    if (mag(eofs.Y()) < SMALL)
    {
        return RDEThermo<EquationOfState>(eofs);
    }
    else
    {
        return RDEThermo<EquationOfState>(eofs);
    }
}

template<class EquationOfState>
inline Foam::RDEThermo<EquationOfState> Foam::operator * (const scalar s, const RDEThermo<EquationOfState>& ct)
{
    return RDEThermo<EquationOfState>(s*static_cast<const EquationOfState&>(ct));
}

template<class EquationOfState>
inline Foam::RDEThermo<EquationOfState> Foam::operator==(const RDEThermo<EquationOfState> & ct1, const RDEThermo<EquationOfState> & ct2)
{
    EquationOfState eofs(static_cast<const EquationOfState&>(ct1) == static_cast<const EquationOfState&>(ct2));
    return RDEThermo<EquationOfState>
    (
        eofs,
        ct2.Y()/eofs.Y()*ct2.Edissociation - ct1.Y()/eofs.Y()*ct1.Edissociation
    );
}

#endif
