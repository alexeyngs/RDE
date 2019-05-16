
#ifndef timeWallProcessorCyclicFvPatchField_H
#define timeWallProcessorCyclicFvPatchField_H

#include "processorCyclicFvPatch.H"
#include "processorFvPatchField.H"
#include "processorCyclicFvPatch.H"
#include "demandDrivenData.H"
#include "transformField.H"
#include "scalarIOList.H"

namespace Foam
{
template<class Type>
class timeWallProcessorCyclicFvPatchField : public processorFvPatchField<Type>
{
    //- Send buffer.
    mutable Field<Type> sendBuffer;
    //- Local reference cast into the processor patch
    const processorCyclicFvPatch & processorPatch;
public:
    //- Runtime type information
    TypeName("timeWallProcessorCyclicVolume");
    //- Construct from patch and internal field
    timeWallProcessorCyclicFvPatchField
    (
        const fvPatch&,
        const DimensionedField<Type, volMesh>&
    );

    //- Construct from patch and internal field and patch field
    timeWallProcessorCyclicFvPatchField
    (
        const fvPatch&,
        const DimensionedField<Type, volMesh>&,
        const Field<Type>&
    );

    //- Construct from patch, internal field and dictionary
    timeWallProcessorCyclicFvPatchField
    (
        const fvPatch&,
        const DimensionedField<Type, volMesh>&,
        const dictionary&
    );

    //- Construct by mapping given timeWallProcessorCyclicFvPatchField onto a
    //  new patch
    timeWallProcessorCyclicFvPatchField
    (
        const timeWallProcessorCyclicFvPatchField<Type>&,
        const fvPatch&,
        const DimensionedField<Type, volMesh>&,
        const fvPatchFieldMapper&
    );

    //- Construct as copy
    timeWallProcessorCyclicFvPatchField(const timeWallProcessorCyclicFvPatchField<Type>&);

    //- Construct and return a clone
    virtual tmp<fvPatchField<Type>> clone() const
    {
        return tmp<fvPatchField<Type>>
        (
            new timeWallProcessorCyclicFvPatchField<Type>(*this)
        );
    }

    //- Construct as copy setting internal field reference
    timeWallProcessorCyclicFvPatchField
    (
        const timeWallProcessorCyclicFvPatchField<Type>&,
        const DimensionedField<Type, volMesh>&
    );

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchField<Type>> clone
    (
        const DimensionedField<Type, volMesh>& iF
    ) const
    {
        return tmp<fvPatchField<Type>>
        (
            new timeWallProcessorCyclicFvPatchField<Type>(*this, iF)
        );
    }
    //- Destructor
    virtual ~timeWallProcessorCyclicFvPatchField();


    // Member functions
    //- Does the patch field perform the transformation
    virtual bool doTransform() const
    {
        return !(processorPatch.parallel() || pTraits<Type>::rank == 0);
    }

    //- Return face transformation tensor
    virtual const tensorField& forwardT() const
    {
        return processorPatch.forwardT();
    }
    // Свои методы
    virtual void initEvaluate(const Pstream::commsTypes commsType) override;
};
} // End namespace Foam




template<class Type>
Foam::timeWallProcessorCyclicFvPatchField<Type>::timeWallProcessorCyclicFvPatchField
(const fvPatch& p, const DimensionedField<Type, volMesh>& iF)
: processorFvPatchField<Type>(p, iF), processorPatch(refCast<const processorCyclicFvPatch>(p))
{}


template<class Type>
Foam::timeWallProcessorCyclicFvPatchField<Type>::timeWallProcessorCyclicFvPatchField
( const fvPatch& p, const DimensionedField<Type, volMesh>& iF, const dictionary& dict)
: processorFvPatchField<Type>(p, iF, dict), processorPatch(refCast<const processorCyclicFvPatch>(p))
{
    if (!isType<processorCyclicFvPatch>(p))
    {
        FatalIOErrorInFunction(dict)   << "\n    patch type '" << p.type() << "' not constraint type '" << typeName << "'"
        << "\n    for patch " << p.name() << " of field " << this->internalField().name() << " in file " << this->internalField().objectPath() << exit(FatalIOError);
    }

    if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
    {
        WarningInFunction << "Scheduled communication with split cyclics not supported." << endl;
    }
}


template<class Type>
Foam::timeWallProcessorCyclicFvPatchField<Type>::timeWallProcessorCyclicFvPatchField
(const fvPatch& p, const DimensionedField<Type, volMesh>& iF, const Field<Type>& f)
: processorFvPatchField<Type>(p, iF, f), processorPatch(refCast<const processorCyclicFvPatch>(p))
{}


template<class Type>
Foam::timeWallProcessorCyclicFvPatchField<Type>::timeWallProcessorCyclicFvPatchField
(const timeWallProcessorCyclicFvPatchField<Type>& ptf, const fvPatch& p, const DimensionedField<Type, volMesh>& iF, const fvPatchFieldMapper& mapper)
: processorFvPatchField<Type>(ptf, p, iF, mapper), processorPatch(refCast<const processorCyclicFvPatch>(p))
{
    if (!isType<processorCyclicFvPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::timeWallProcessorCyclicFvPatchField<Type>::timeWallProcessorCyclicFvPatchField
(const timeWallProcessorCyclicFvPatchField<Type>& ptf)
: processorFvPatchField<Type>(ptf), processorPatch(refCast<const processorCyclicFvPatch>(ptf.patch()))
{}


template<class Type>
Foam::timeWallProcessorCyclicFvPatchField<Type>::timeWallProcessorCyclicFvPatchField
(const timeWallProcessorCyclicFvPatchField<Type>& ptf, const DimensionedField<Type, volMesh>& iF)
: processorFvPatchField<Type>(ptf, iF), processorPatch(refCast<const processorCyclicFvPatch>(ptf.patch()))
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class Type>
Foam::timeWallProcessorCyclicFvPatchField<Type>::~timeWallProcessorCyclicFvPatchField()
{}

//--------------------------------------------------------------------------------
template<class Type>
void Foam::timeWallProcessorCyclicFvPatchField<Type>::initEvaluate(const Pstream::commsTypes commsType)
{
    // Получаем список скаляргых параметров, которые хранятся в главное программе
    const Foam::scalarIOList & scalarParameters = this->db().template lookupObject<scalarIOList>("scalarParameters");
    // Из этого списка берем нужное
    const Foam::scalar & timeOnePeriod = scalarParameters[0];
    const Foam::scalar & time = this->db().time().value();   // Время
    if(time > timeOnePeriod * 0.9) // Периодическое условие
    {
        processorFvPatchField<Type>::initEvaluate(commsType);
    }else
    {
        this->patchInternalField(sendBuffer);
        this->operator == (sendBuffer);
    }
}
#endif
