
#ifndef RDEspecie_H
#define RDEspecie_H
#include "word.H"
#include "scalar.H"
#include "dictionary.H"
#include "thermodynamicConstants.H"
using namespace Foam::constant::thermodynamic;

namespace Foam
{
// Forward declaration of friend functions and operators
class RDEspecie;
inline RDEspecie operator + (const RDEspecie&, const RDEspecie&);
inline RDEspecie operator * (const scalar, const RDEspecie&);
inline RDEspecie operator == (const RDEspecie&, const RDEspecie&);
Ostream & operator << (Ostream&, const RDEspecie&);

class RDEspecie
{
public:

        //- Name of specie
        word name_;

        //- Number of moles of this component in the mixture
        scalar moles_;

        // - Количество атомов в смеси
        scalar atoms_;

        //- Количество атомов, на которых будет разлагаться смесь
        scalar decomposition_;

        //- Количество атомов кислорода, которое необходимо для полного сжигания компонента
        scalar needO_;

        //- Молекулярная масса [kg/kmol]
        mutable scalar molWeight_;

        //- Энергия диссоциации компонента смеси        
        mutable scalar Edissociation_;

        //- Энергия активации [Дж/моль]
        mutable scalar Einduction_;

        //- Константа активации в моль*с/м3
        mutable scalar Kinduction_;
public:
    //- Runtime type information
    ClassName("RDEspecie");
    // Constructors
    //- Construct from components without name
    inline RDEspecie(const scalar moles, const scalar atoms, const scalar decomposition, const scalar molWeight, const scalar needO, const scalar Edissociation, const scalar Einduction, const scalar Kinduction);

    //- Construct from components with name
    inline RDEspecie(const word & name, const scalar moles, const scalar atoms, const scalar decomposition, const scalar molWeight, const scalar needO, const scalar Edissociation, const scalar Einduction, const scalar Kinduction);

    //- Construct as named copy
    inline RDEspecie(const word & name, const RDEspecie&);

    //- Construct from dictionary
    RDEspecie(const dictionary & dict);
// Member Functions
    //- Name
    inline const word & name() const;

    //- Молекулярная масса [kg/kmol]
    inline scalar W() const;

    //- Количество атомов кислорода, которое необходимо для полного сжигания компонента
    inline scalar needO() const;

    //- Количество молей
    inline scalar moles() const;

    //- Количество атомов в смеси
    inline scalar atoms() const;

    //- Количество атомов, на которых будет разлагаться смесь
    inline scalar decomposition() const;

    //- Энергия диссоциации компонента смеси        
    inline scalar Edissociation() const;

    //- Энергия активации [Дж/моль]
    inline scalar Einduction() const;

    //- Константа активации в моль*с/м3
    inline scalar Kinduction() const;

    //- No of moles of this RDEspecies in mixture
    inline scalar Y() const;

    //- Gas constant [J/(kg K)]
    inline scalar R() const;

//- Write to Ostream
    void write(Ostream& os) const;
    // Member operators
    inline void operator=(const RDEspecie&);
    inline void operator+=(const RDEspecie&);
    inline void operator*=(const scalar);

    // Friend operators
    inline friend RDEspecie operator+(const RDEspecie&, const RDEspecie&);
    inline friend RDEspecie operator*(const scalar, const RDEspecie&);
    inline friend RDEspecie operator==(const RDEspecie&, const RDEspecie&);
    // Ostream Operator
    friend Ostream& operator<<(Ostream&, const RDEspecie&);
};

// Конструкторы
inline RDEspecie::RDEspecie(const word & name, const scalar moles, const scalar atoms, const scalar decomposition, const scalar molWeight,
    const scalar needO, const scalar Edissociation, const scalar Einduction, const scalar Kinduction):
    name_(name), moles_(moles), atoms_(atoms), decomposition_(decomposition), molWeight_(molWeight), needO_(needO),
    Edissociation_(Edissociation), Einduction_(Einduction), Kinduction_(Kinduction)
{}

inline RDEspecie::RDEspecie(const scalar moles, const scalar atoms, const scalar decomposition, const scalar molWeight,
    const scalar needO, const scalar Edissociation, const scalar Einduction, const scalar Kinduction):
    moles_(moles), atoms_(atoms), decomposition_(decomposition), molWeight_(molWeight), needO_(needO),
    Edissociation_(Edissociation), Einduction_(Einduction), Kinduction_(Kinduction)
{}

inline RDEspecie::RDEspecie(const word & name, const RDEspecie & st):
    name_(name), moles_(st.moles_), molWeight_(st.molWeight_), needO_(st.needO_), atoms_(st.atoms_),
    Edissociation_(st.Edissociation_), Einduction_(st.Einduction_), Kinduction_(st.Kinduction_),
    decomposition_(st.decomposition_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const word& RDEspecie::name() const
{
    return name_;
}

inline scalar RDEspecie::W() const
{
    return molWeight_;
}

inline scalar RDEspecie::moles() const
{
    return moles_;
}

inline scalar RDEspecie::atoms() const
{
    return atoms_;
}

inline scalar RDEspecie::decomposition() const
{
    return decomposition_;
}

inline scalar RDEspecie::Edissociation() const
{
    return Edissociation_;
}

inline scalar RDEspecie::Einduction() const
{
    return Einduction_;
}

inline scalar RDEspecie::Kinduction() const
{
    return Kinduction_;
}

inline scalar RDEspecie::needO() const
{
    return needO_;
}

inline scalar RDEspecie::Y() const
{
    return 1.0;
}


inline scalar RDEspecie::R() const
{
    return RR/molWeight_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void RDEspecie::operator=(const RDEspecie & st)
{
    //name_ = st.name_;
    moles_ = st.moles_;
    atoms_ = st.atoms_;
    decomposition_ = st.decomposition_;
    molWeight_ = st.molWeight_;
    needO_ = st.needO_;
    Edissociation_ = st.Edissociation_;
    Einduction_ = st.Einduction_;
    Kinduction_ = st.Kinduction_;
}


inline void RDEspecie::operator += (const RDEspecie & st)
{
    const scalar summoles = moles_ + st.moles_;
    // st.moles_ - сумма старых молярных масс
    // summoles - сумма новых молярных масс
    if (mag(summoles) > SMALL)
    {
        atoms_ = (moles_*atoms_ + st.moles_*st.atoms_) / summoles;
        decomposition_ = (moles_*decomposition_ + st.moles_*st.decomposition_) / summoles;
        molWeight_ = (moles_*molWeight_ + st.moles_*st.molWeight_) / summoles;
        needO_ = (moles_*needO_ + st.moles_*st.needO_) / summoles;
        Edissociation_ = (moles_*Edissociation_ + st.moles_*st.Edissociation_) / summoles;
        Einduction_ = (moles_*Einduction_ + st.moles_*st.Einduction_) / summoles;
        Kinduction_ = (moles_*Kinduction_ + st.moles_*st.Kinduction_) / summoles;
    }
    moles_ = summoles;
}

inline void RDEspecie::operator *= (const scalar s)
{
    moles_ *= s;
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline RDEspecie operator + (const RDEspecie& st1, const RDEspecie& st2)
{
    const scalar summoles = max(st1.moles_ + st2.moles_, SMALL);

    if (mag(summoles) > SMALL)
    {
        return RDEspecie
        (
            summoles,
            (st1.moles_*st1.atoms_ + st2.moles_*st2.atoms_) / summoles,
            (st1.moles_*st1.decomposition_ + st2.moles_*st2.decomposition_) / summoles,
            (st1.moles_*st1.molWeight_ + st2.moles_*st2.molWeight_) / summoles,
            (st1.moles_*st1.needO_ + st2.moles_*st2.needO_) / summoles,
            (st1.moles_*st1.Edissociation_ + st2.moles_*st2.Edissociation_) / summoles,
            (st1.moles_*st1.Einduction_ + st2.moles_*st2.Einduction_) / summoles,
            (st1.moles_*st1.Kinduction_ + st2.moles_*st2.Kinduction_) / summoles
        );
    }
    else
    {
        return st1;
    }
}

// Зачем-то сюда заходит...
inline RDEspecie operator * (const scalar s, const RDEspecie & st)
{
    return RDEspecie
    (
        st.moles_,
        st.atoms_,
        st.decomposition_,
        st.molWeight_,
        st.needO_,
        st.Edissociation_,
        st.Einduction_,
        st.Kinduction_
    );
}


inline RDEspecie operator==(const RDEspecie & st1, const RDEspecie & st2)
{
    scalar diff_moles = st2.moles_ - st1.moles_;
    if (mag(diff_moles) < SMALL)
    {
        diff_moles = SMALL;
    }

    const scalar diff_atoms = st2.moles_*st2.atoms_ - st1.moles_*st1.atoms_;
    const scalar diff_decomposition = st2.moles_*st2.decomposition_ - st1.moles_*st1.decomposition_;
    const scalar diff_molWeight = st2.moles_*st2.molWeight_ - st1.moles_*st1.molWeight_;
    const scalar diff_NeedO = st2.moles_*st2.needO_ - st1.moles_*st1.needO_;
    const scalar diff_Edissociation = st2.moles_*st2.Edissociation_ - st1.moles_*st1.Edissociation_;
    const scalar diff_Einduction = st2.moles_*st2.Einduction_ - st1.moles_*st1.Einduction_;
    const scalar diff_Kinduction = st2.moles_*st2.Kinduction_ - st1.moles_*st1.Kinduction_;
    
    scalar atoms = 0.0;
    scalar decomposition = 0.0;
    scalar molWeight = 0.0;
    scalar needO = 0.0;
    scalar Edissociation = 0.0;
    scalar Einduction = 0.0;
    scalar Kinduction = 0.0;

    if (mag(diff_atoms) > SMALL) atoms = diff_atoms / diff_moles;
    if (mag(diff_decomposition) > SMALL) decomposition = diff_decomposition / diff_moles;
    if (mag(diff_molWeight) > SMALL) molWeight = diff_molWeight / diff_moles;
    if (mag(diff_NeedO) > SMALL) needO = diff_NeedO / diff_moles;
    if (mag(diff_Edissociation) > SMALL) Edissociation = diff_Edissociation / diff_moles;
    if (mag(diff_Einduction) > SMALL) Einduction = diff_Einduction / diff_moles;
    if (mag(diff_Kinduction) > SMALL) Kinduction = diff_Kinduction / diff_moles;

    return RDEspecie(diff_moles, atoms, decomposition, molWeight, needO, Edissociation, Einduction, Kinduction);
}

} // End namespace Foam
#endif
