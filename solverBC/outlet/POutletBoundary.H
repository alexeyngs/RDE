
#ifndef POutletBoundary_H
#define POutletBoundary_H
#include "fvPatchField.H"
#include "fixedValueFvPatchFields.H"
#include "boundary.H"
namespace Foam
{
class POutletBoundary : public fixedValueFvPatchField<scalar>
{
public:
    //- Information about the SHA1 of the code itself
    static const char* const SHA1sum;
    //- Runtime type information
    TypeName("POutletBoundary");
    // Constructors
    //- Construct from patch and internal field
    POutletBoundary
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct from patch, internal field and dictionary
    POutletBoundary(const fvPatch&,const DimensionedField<scalar, volMesh>&,const dictionary&);
    //- Construct by mapping a copy onto a new patch
    POutletBoundary(const POutletBoundary&,const fvPatch&,const DimensionedField<scalar, volMesh>&,const fvPatchFieldMapper&);
    //- Construct as copy
    POutletBoundary(const POutletBoundary&);

    //- Construct and return a clone
    virtual tmp<fvPatchScalarField > clone() const
    {
        return tmp<fvPatchScalarField> (new POutletBoundary(*this));
    }

    //- Construct as copy setting internal field reference
    POutletBoundary(const POutletBoundary&,const DimensionedField<scalar, volMesh>&);

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchScalarField > clone(const DimensionedField<scalar, volMesh>& iF) const
    {
        return tmp<fvPatchScalarField >
        (
            new POutletBoundary(*this, iF)
        );
    }
    //- Destructor
    virtual ~POutletBoundary();
    // Member functions
    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();
};
} // End namespace Foam
#endif
