#ifndef timeWallProcessorCyclicFvsPatchField_H
#define timeWallProcessorCyclicFvsPatchField_H

#include "coupledFvsPatchField.H"
#include "processorCyclicFvPatch.H"
namespace Foam
{
template<class Type>
class timeWallProcessorCyclicFvsPatchField : public coupledFvsPatchField<Type>
{
    //- Local reference cast into the processor patch
    const processorCyclicFvPatch & processorPatch;
public:

    //- Runtime type information
    TypeName("timeWallProcessorCyclicSurface");
    //- Construct from patch and internal field
    timeWallProcessorCyclicFvsPatchField
    (
        const fvPatch&,
        const DimensionedField<Type, surfaceMesh>&
    );

    //- Construct from patch and internal field and patch field
    timeWallProcessorCyclicFvsPatchField
    (
        const fvPatch&,
        const DimensionedField<Type, surfaceMesh>&,
        const Field<Type>&
    );

    //- Construct from patch, internal field and dictionary
    timeWallProcessorCyclicFvsPatchField
    (
        const fvPatch&,
        const DimensionedField<Type, surfaceMesh>&,
        const dictionary&
    );

    //- Construct by mapping given timeWallProcessorCyclicFvsPatchField onto a
    //  new patch
    timeWallProcessorCyclicFvsPatchField
    (
        const timeWallProcessorCyclicFvsPatchField<Type>&,
        const fvPatch&,
        const DimensionedField<Type, surfaceMesh>&,
        const fvPatchFieldMapper&
    );

    //- Construct as copy
    timeWallProcessorCyclicFvsPatchField(const timeWallProcessorCyclicFvsPatchField<Type>&);

    //- Construct and return a clone
    virtual tmp<fvsPatchField<Type>> clone() const
    {
        return tmp<fvsPatchField<Type>>
        (
            new timeWallProcessorCyclicFvsPatchField<Type>(*this)
        );
    }

    //- Construct as copy setting internal field reference
    timeWallProcessorCyclicFvsPatchField
    ( const timeWallProcessorCyclicFvsPatchField<Type> &, const DimensionedField<Type, surfaceMesh> &);

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvsPatchField<Type>> clone(const DimensionedField<Type, surfaceMesh>& iF) const
    {
        return tmp<fvsPatchField<Type>>
        (
            new timeWallProcessorCyclicFvsPatchField<Type>(*this, iF)
        );
    }


    //- Destructor
    virtual ~timeWallProcessorCyclicFvsPatchField();
    // Member functions

    //- Return true if running parallel
    virtual bool coupled() const
    {
        if (Pstream::parRun())
        {
            return true;
        }
        else
        {
            return false;
        }
    }
};
} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Type>
Foam::timeWallProcessorCyclicFvsPatchField<Type>::timeWallProcessorCyclicFvsPatchField
(const fvPatch & p, const DimensionedField<Type, surfaceMesh>& iF)
: coupledFvsPatchField<Type>(p, iF), processorPatch(refCast<const processorCyclicFvPatch>(p))
{
}


template<class Type>
Foam::timeWallProcessorCyclicFvsPatchField<Type>::timeWallProcessorCyclicFvsPatchField
(const fvPatch& p, const DimensionedField<Type, surfaceMesh>& iF, const Field<Type>& f)
: coupledFvsPatchField<Type>(p, iF, f), processorPatch(refCast<const processorCyclicFvPatch>(p))
{}


template<class Type>
Foam::timeWallProcessorCyclicFvsPatchField<Type>::timeWallProcessorCyclicFvsPatchField
(const fvPatch & p, const DimensionedField<Type, surfaceMesh>& iF, const dictionary& dict)
: coupledFvsPatchField<Type>(p, iF, dict), processorPatch(refCast<const processorCyclicFvPatch>(p))
{
    if (!isType<processorCyclicFvPatch>(p))
    {
        FatalIOErrorInFunction(dict)   << "patch " << this->patch().index() << " not processor type. "
        << "Patch type = " << p.type() << exit(FatalIOError);
    }
}


template<class Type>
Foam::timeWallProcessorCyclicFvsPatchField<Type>::timeWallProcessorCyclicFvsPatchField
(const timeWallProcessorCyclicFvsPatchField<Type>& ptf, const fvPatch& p, const DimensionedField<Type, surfaceMesh>& iF, const fvPatchFieldMapper& mapper)
: coupledFvsPatchField<Type>(ptf, p, iF, mapper), processorPatch(refCast<const processorCyclicFvPatch>(p))
{
    if (!isType<processorCyclicFvPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template<class Type>
Foam::timeWallProcessorCyclicFvsPatchField<Type>::timeWallProcessorCyclicFvsPatchField
(const timeWallProcessorCyclicFvsPatchField<Type>& ptf)
: coupledFvsPatchField<Type>(ptf), processorPatch(refCast<const processorCyclicFvPatch>(ptf.patch()))
{}


template<class Type>
Foam::timeWallProcessorCyclicFvsPatchField<Type>::timeWallProcessorCyclicFvsPatchField
(const timeWallProcessorCyclicFvsPatchField<Type>& ptf, const DimensionedField<Type, surfaceMesh>& iF)
: coupledFvsPatchField<Type>(ptf, iF), processorPatch(refCast<const processorCyclicFvPatch>(ptf.patch()))
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class Type>
Foam::timeWallProcessorCyclicFvsPatchField<Type>::~timeWallProcessorCyclicFvsPatchField()
{}

#endif
