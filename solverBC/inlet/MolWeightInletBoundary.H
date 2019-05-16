
#ifndef MolWeightInletBoundary_H
#define MolWeightInletBoundary_H
#include "fvPatchField.H"
#include "fixedValueFvPatchFields.H"
#include "boundary.H"
namespace Foam
{
    
class MolWeightInletBoundary : public fixedValueFvPatchField<scalar>
{
public:
    //- Information about the SHA1 of the code itself
    static const char* const SHA1sum;
    //- Runtime type information
    TypeName("MolWeightInletBoundary");
    // Constructors
    //- Construct from patch and internal field
    MolWeightInletBoundary
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct from patch, internal field and dictionary
    MolWeightInletBoundary(const fvPatch&,const DimensionedField<scalar, volMesh>&,const dictionary&);
    //- Construct by mapping a copy onto a new patch
    MolWeightInletBoundary(const MolWeightInletBoundary&,const fvPatch&,const DimensionedField<scalar, volMesh>&,const fvPatchFieldMapper&);
    //- Construct as copy
    MolWeightInletBoundary(const MolWeightInletBoundary&);

    //- Construct and return a clone
    virtual tmp<fvPatchScalarField > clone() const
    {
        return tmp<fvPatchScalarField> (new MolWeightInletBoundary(*this));
    }

    //- Construct as copy setting internal field reference
    MolWeightInletBoundary(const MolWeightInletBoundary&,const DimensionedField<scalar, volMesh>&);

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchScalarField > clone(const DimensionedField<scalar, volMesh>& iF) const
    {
        return tmp<fvPatchScalarField >
        (
            new MolWeightInletBoundary(*this, iF)
        );
    }
    //- Destructor
    virtual ~MolWeightInletBoundary();
    // Member functions
    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();
};
} // End namespace Foam
#endif
