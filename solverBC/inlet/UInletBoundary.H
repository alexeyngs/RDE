#ifndef UInletBoundary_H
#define UInletBoundary_H
#include "fvPatchField.H"
#include "fixedValueFvPatchFields.H"
#include "boundary.H"
namespace Foam
{
class UInletBoundary : public fixedValueFvPatchField<vector>
{
public:
    //- Information about the SHA1 of the code itself
    static const char* const SHA1sum;

    //- Runtime type information
    TypeName("UInletBoundary");

    // Constructors
    //- Construct from patch and internal field
    UInletBoundary(const fvPatch&,const DimensionedField<vector, volMesh>&);

    //- Construct from patch, internal field and dictionary
    UInletBoundary(const fvPatch&,const DimensionedField<vector, volMesh>&,const dictionary&);

    //- Construct by mapping a copy onto a new patch
    UInletBoundary(const UInletBoundary&,const fvPatch&,const DimensionedField<vector, volMesh>&,const fvPatchFieldMapper&);

    //- Construct as copy
    UInletBoundary(const UInletBoundary&);

    //- Construct and return a clone
    virtual tmp<fvPatchVectorField > clone() const
    {
        return tmp<fvPatchVectorField >
        (
            new UInletBoundary(*this)
        );
    }

    //- Construct as copy setting internal field reference
    UInletBoundary(const UInletBoundary&,const DimensionedField<vector, volMesh>&);

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchVectorField > clone(const DimensionedField<vector, volMesh>& iF) const
    {
        return tmp<fvPatchVectorField >
        (
            new UInletBoundary(*this, iF)
        );
    }
    //- Destructor
    virtual ~UInletBoundary();
    // Member functions
    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();
};
} // End namespace Foam
#endif