#include "MolWeightInletBoundary.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
#include "thermodynamicConstants.H"
#include "basicSpecieMixture.H"
#include "rhoReactionThermo.H"
#include "../solver/thermo/rhoRDEThermo.H"

namespace Foam
{
// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //
//{{{ begin localCode
//}}} end localCode
// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //
extern "C"
{
    // unique function name that can be checked if the correct library version has been loaded
    void MolWeightInletBoundary_1_0(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchScalarField,
    MolWeightInletBoundary
);
const char * const MolWeightInletBoundary::SHA1sum = "66d4e0b38a18f2b7674332f6d021721a92f0e80e";

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
MolWeightInletBoundary::MolWeightInletBoundary(const fvPatch & p, const DimensionedField<scalar, volMesh> & iF):
fixedValueFvPatchField<scalar>(p, iF)
{
    if (false)
    {
        Info<<"construct MolWeightInletBoundary sha1: 66d4e0b38a18f2b7674332f6d021721a92f0e80e from patch/DimensionedField" << endl;
    }
}
MolWeightInletBoundary::MolWeightInletBoundary(const MolWeightInletBoundary& ptf,const fvPatch& p,const DimensionedField<scalar, volMesh>& iF,const fvPatchFieldMapper& mapper):
fixedValueFvPatchField<scalar>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct MolWeightInletBoundary sha1: 66d4e0b38a18f2b7674332f6d021721a92f0e80e from patch/DimensionedField/mapper\n";
    }
}
MolWeightInletBoundary::MolWeightInletBoundary(const fvPatch& p,const DimensionedField<scalar, volMesh>& iF,const dictionary& dict):
fixedValueFvPatchField<scalar>(p, iF, dict)
{
    if (false)
    {
        Info <<"construct MolWeightInletBoundary sha1: 66d4e0b38a18f2b7674332f6d021721a92f0e80e from patch/dictionary\n";
    }
}

MolWeightInletBoundary::MolWeightInletBoundary(const MolWeightInletBoundary& ptf):
fixedValueFvPatchField<scalar>(ptf)
{
    if (false)
    {
        Info <<"construct MolWeightInletBoundary sha1: 66d4e0b38a18f2b7674332f6d021721a92f0e80e as copy\n";
    }
}

MolWeightInletBoundary::MolWeightInletBoundary(const MolWeightInletBoundary& ptf,const DimensionedField<scalar, volMesh>& iF):
fixedValueFvPatchField<scalar>(ptf, iF)
{
    if (false)
    {
        Info <<"construct MolWeightInletBoundary sha1: 66d4e0b38a18f2b7674332f6d021721a92f0e80e as copy/DimensionedField\n";
    }
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
MolWeightInletBoundary::~MolWeightInletBoundary()
{
    if (false)
    {
        Info <<"destroy MolWeightInletBoundary sha1: 66d4e0b38a18f2b7674332f6d021721a92f0e80e\n";
    }
}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void MolWeightInletBoundary::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    if (false)
    {
        Info<<"updateCoeffs MolWeightInletBoundary sha1: 66d4e0b38a18f2b7674332f6d021721a92f0e80e\n";
    }
// begin code
    const scalar & time = this->db().time().value();   // Время
    const Foam::scalarIOList & scalarParameters = this->db().template lookupObject<scalarIOList>("scalarParameters");
    const Foam::vectorField & Cf = this->patch().Cf();
    // Получим текущий патч - границу
    const label ThisPatch = patch().index();
    scalarField & result = *this;
    // на граничной ячейки
    tmp<scalarField> Field(this->patchInternalField());
    tmp<scalarField> rho = this->db().template lookupObject<volScalarField>("rho").boundaryField()[ThisPatch].patchInternalField();
    tmp<scalarField> p = this->db().template lookupObject<volScalarField>("p").boundaryField()[ThisPatch].patchInternalField();
    tmp<vectorField> U = this->db().template lookupObject<volVectorField>("U").boundaryField()[ThisPatch].patchInternalField();
    tmp<scalarField> gamma = this->db().template lookupObject<volScalarField>("gamma").boundaryField()[ThisPatch].patchInternalField();
    //-----------------------------------------------------
    const rhoReactionThermo & thermo = this->db().lookupObject<rhoReactionThermo>("thermophysicalProperties");
    rhoRDEThermo & Thermo = (rhoRDEThermo &)thermo;
    forAll(Cf, facei)
    {
        result[facei] = Thermo.patchMolWeight(ThisPatch, facei);
    }
// end code
    this->fixedValueFvPatchField<scalar>::updateCoeffs();
}
} // End namespace Foam
