
#ifndef MolWeightOutletBoundary_H
#define MolWeightOutletBoundary_H
#include "fvPatchField.H"
#include "fixedValueFvPatchFields.H"
#include "boundary.H"
namespace Foam
{
    
class MolWeightOutletBoundary : public fixedValueFvPatchField<scalar>
{
public:
    //- Information about the SHA1 of the code itself
    static const char* const SHA1sum;
    //- Runtime type information
    TypeName("MolWeightOutletBoundary");
    // Constructors
    //- Construct from patch and internal field
    MolWeightOutletBoundary
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct from patch, internal field and dictionary
    MolWeightOutletBoundary(const fvPatch&,const DimensionedField<scalar, volMesh>&,const dictionary&);
    //- Construct by mapping a copy onto a new patch
    MolWeightOutletBoundary(const MolWeightOutletBoundary&,const fvPatch&,const DimensionedField<scalar, volMesh>&,const fvPatchFieldMapper&);
    //- Construct as copy
    MolWeightOutletBoundary(const MolWeightOutletBoundary&);

    //- Construct and return a clone
    virtual tmp<fvPatchScalarField > clone() const
    {
        return tmp<fvPatchScalarField> (new MolWeightOutletBoundary(*this));
    }

    //- Construct as copy setting internal field reference
    MolWeightOutletBoundary(const MolWeightOutletBoundary&,const DimensionedField<scalar, volMesh>&);

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchScalarField > clone(const DimensionedField<scalar, volMesh>& iF) const
    {
        return tmp<fvPatchScalarField >
        (
            new MolWeightOutletBoundary(*this, iF)
        );
    }
    //- Destructor
    virtual ~MolWeightOutletBoundary();
    // Member functions
    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();
};
} // End namespace Foam
#endif
