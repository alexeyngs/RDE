#ifndef UOutletBoundary_H
#define UOutletBoundary_H
#include "fvPatchField.H"
#include "fixedValueFvPatchFields.H"
#include "boundary.H"
namespace Foam
{
class UOutletBoundary : public fixedValueFvPatchField<vector>
{
public:
    //- Information about the SHA1 of the code itself
    static const char* const SHA1sum;

    //- Runtime type information
    TypeName("UOutletBoundary");

    // Constructors
    //- Construct from patch and internal field
    UOutletBoundary(const fvPatch&,const DimensionedField<vector, volMesh>&);

    //- Construct from patch, internal field and dictionary
    UOutletBoundary(const fvPatch&,const DimensionedField<vector, volMesh>&,const dictionary&);

    //- Construct by mapping a copy onto a new patch
    UOutletBoundary(const UOutletBoundary&,const fvPatch&,const DimensionedField<vector, volMesh>&,const fvPatchFieldMapper&);

    //- Construct as copy
    UOutletBoundary(const UOutletBoundary&);

    //- Construct and return a clone
    virtual tmp<fvPatchVectorField > clone() const
    {
        return tmp<fvPatchVectorField >
        (
            new UOutletBoundary(*this)
        );
    }

    //- Construct as copy setting internal field reference
    UOutletBoundary(const UOutletBoundary&,const DimensionedField<vector, volMesh>&);

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchVectorField > clone(const DimensionedField<vector, volMesh>& iF) const
    {
        return tmp<fvPatchVectorField >
        (
            new UOutletBoundary(*this, iF)
        );
    }
    //- Destructor
    virtual ~UOutletBoundary();
    // Member functions
    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();
};
} // End namespace Foam
#endif