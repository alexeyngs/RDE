#ifndef timeWallCyclicFvPatchField_H
#define timeWallCyclicFvPatchField_H

#include "cyclicFvPatchField.H"
#include "scalarIOList.H"
namespace Foam
{

template<class Type>
class timeWallCyclicFvPatchField : public cyclicFvPatchField<Type>
{
public:
    TypeName("timeWallCyclic");
    // Constructors

    //- Construct from patch and internal field
    timeWallCyclicFvPatchField(const fvPatch&,const DimensionedField<Type, volMesh>&);

    //- Construct from patch, internal field and dictionary
    timeWallCyclicFvPatchField(const fvPatch&, const DimensionedField<Type, volMesh>&, const dictionary&);

    //- Construct by mapping given timeWallCyclicFvPatchField onto a new patch
    timeWallCyclicFvPatchField(const timeWallCyclicFvPatchField<Type>&, const fvPatch&, const DimensionedField<Type, volMesh>&, const fvPatchFieldMapper&);

    //- Construct as copy
    timeWallCyclicFvPatchField(const timeWallCyclicFvPatchField<Type>&);

    //- Construct and return a clone
    virtual tmp<fvPatchField<Type>> clone() const
    {
        return tmp<fvPatchField<Type>> (new timeWallCyclicFvPatchField<Type>(*this));
    }

    //- Construct as copy setting internal field reference
    timeWallCyclicFvPatchField(const timeWallCyclicFvPatchField<Type>&,const DimensionedField<Type, volMesh>&);

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchField<Type>> clone(const DimensionedField<Type, volMesh>& iF) const
    {
        return tmp<fvPatchField<Type>>
        (
            new timeWallCyclicFvPatchField<Type>(*this, iF)
        );
    }
    // Member functions
    //- Return neighbour coupled internal cell data
    tmp<Field<Type>> patchNeighbourField() const;
    virtual void evaluate(const Pstream::commsTypes stream) override;
    void write(Ostream& os) const;
};
} // End namespace Foam
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<class Type>
Foam::timeWallCyclicFvPatchField<Type>::timeWallCyclicFvPatchField(const fvPatch& p,const DimensionedField<Type, volMesh>& iF):cyclicFvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::timeWallCyclicFvPatchField<Type>::timeWallCyclicFvPatchField
(
    const timeWallCyclicFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
):cyclicFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::timeWallCyclicFvPatchField<Type>::timeWallCyclicFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
):cyclicFvPatchField<Type>(p, iF, dict)
{}


template<class Type>
Foam::timeWallCyclicFvPatchField<Type>::timeWallCyclicFvPatchField
(
    const timeWallCyclicFvPatchField<Type>& ptf
):cyclicFvPatchField<Type>(ptf)
{}


template<class Type>
Foam::timeWallCyclicFvPatchField<Type>::timeWallCyclicFvPatchField
(
    const timeWallCyclicFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
):cyclicFvPatchField<Type>(ptf, iF)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class TScalar> inline TScalar ABS(TScalar X)
{
	return X > 0 ? X : -X;
};

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::timeWallCyclicFvPatchField<Type>::patchNeighbourField() const
{
    const Foam::Field<Type>& iField = this->primitiveField();
    const Foam::labelUList & neighbourFaceCells = cyclicFvPatchField<Type>::cyclicPatch().cyclicPatch().neighbPatch().faceCells();
    const Foam::labelUList & ownerFaceCells = cyclicFvPatchField<Type>::cyclicPatch().cyclicPatch().faceCells();
    // Получаем список скаляргых параметров, которые хранятся в главное программе
    const Foam::scalarIOList & scalarParameters = this->db().template lookupObject<scalarIOList>("scalarParameters");
    // Из этого списка берем нужное
    const Foam::scalar & timeOnePeriod = scalarParameters[0];
    const Foam::scalar & time = this->db().time().value();   // Время

    Foam::tmp<Field<Type>> tpnf(this->patchInternalField());
    Foam::Field<Type> & pnf = tpnf.ref();

    if(time > timeOnePeriod * 0.9) // Периодическое условие
    {
        if (cyclicFvPatchField<Type>::doTransform())
        {
            forAll(pnf, facei)
            {
                pnf[facei] = transform(cyclicFvPatchField<Type>::forwardT()[0], iField[neighbourFaceCells[facei]]);
            }
        }
        else
        {
            forAll(pnf, facei)
            {
                pnf[facei] = iField[neighbourFaceCells[facei]];
            }
        }
    } // Иначе - нулевой градиент
    return tpnf;
}

template<class Type>
void Foam::timeWallCyclicFvPatchField<Type>::evaluate(const Pstream::commsTypes stream)
{
    // Получаем список скаляргых параметров, которые хранятся в главное программе
    const Foam::scalarIOList & scalarParameters = this->db().template lookupObject<scalarIOList>("scalarParameters");
    // Из этого списка берем нужное
    const Foam::scalar & timeOnePeriod = scalarParameters[0];
    const Foam::scalar & time = this->db().time().value();   // Время
    if(time > timeOnePeriod * 0.9) // Периодическое условие
    {
        cyclicFvPatchField<Type>::evaluate(stream);
    }else
    {
        //Info << "#" << endl;
        Field<Type>::operator = (this->patchInternalField());
        fvPatchField<Type>::evaluate();
    }
}

template<class Type>
void Foam::timeWallCyclicFvPatchField<Type>::write(Ostream & os) const
{
    cyclicFvPatchField<Type>::write(os);
    this->writeEntry("value", os);
}
#endif
