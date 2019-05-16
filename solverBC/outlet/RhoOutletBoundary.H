
#ifndef RhoOutletBoundary_H
#define RhoOutletBoundary_H
#include "fvPatchField.H"
#include "fixedValueFvPatchFields.H"
#include "boundary.H"
namespace Foam
{
class RhoOutletBoundary : public fixedValueFvPatchField<scalar>
{
public:
    //- Information about the SHA1 of the code itself
    static const char* const SHA1sum;
    //- Runtime type information
    TypeName("RhoOutletBoundary");
    // Constructors
    //- Construct from patch and internal field
    RhoOutletBoundary
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct from patch, internal field and dictionary
    RhoOutletBoundary(const fvPatch&,const DimensionedField<scalar, volMesh>&,const dictionary&);
    //- Construct by mapping a copy onto a new patch
    RhoOutletBoundary(const RhoOutletBoundary&,const fvPatch&,const DimensionedField<scalar, volMesh>&,const fvPatchFieldMapper&);
    //- Construct as copy
    RhoOutletBoundary(const RhoOutletBoundary&);

    //- Construct and return a clone
    virtual tmp<fvPatchScalarField > clone() const
    {
        return tmp<fvPatchScalarField> (new RhoOutletBoundary(*this));
    }

    //- Construct as copy setting internal field reference
    RhoOutletBoundary(const RhoOutletBoundary&,const DimensionedField<scalar, volMesh>&);

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchScalarField > clone(const DimensionedField<scalar, volMesh>& iF) const
    {
        return tmp<fvPatchScalarField >
        (
            new RhoOutletBoundary(*this, iF)
        );
    }
    //- Destructor
    virtual ~RhoOutletBoundary();
    // Member functions
    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();
};
} // End namespace Foam
#endif
