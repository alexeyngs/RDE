#include "UOutletBoundary.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
#include "thermodynamicConstants.H"
#include "rhoReactionThermo.H"
#include "../solver/thermo/rhoRDEThermo.H"
namespace Foam
{
extern "C"
{
    // unique function name that can be checked if the correct library version has been loaded
    void UOutletBoundary_1_0(bool load)
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

makeRemovablePatchTypeField(fvPatchVectorField,UOutletBoundary);
const char* const UOutletBoundary::SHA1sum = "2cb7a0c4345a61b88d329a26e832e625c53cdd92";

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
UOutletBoundary::UOutletBoundary(const fvPatch& p,const DimensionedField<vector, volMesh>& iF):fixedValueFvPatchField<vector>(p, iF)
{
    if (false)
    {
        Info<<"construct ramU sha1: 2cb7a0c4345a61b88d329a26e832e625c53cdd92"
            " from patch/DimensionedField\n";
    }
}

UOutletBoundary::UOutletBoundary(const UOutletBoundary& ptf,const fvPatch& p,const DimensionedField<vector, volMesh>& iF,const fvPatchFieldMapper& mapper):fixedValueFvPatchField<vector>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct from patch/DimensionedField/mapper\n";
    }
}

UOutletBoundary::UOutletBoundary(const fvPatch& p,const DimensionedField<vector, volMesh>& iF,const dictionary& dict):fixedValueFvPatchField<vector>(p, iF, dict)
{
    if (false)
    {
        Info<<"construct ramU sha1: 2cb7a0c4345a61b88d329a26e832e625c53cdd92 from patch/dictionary\n";
    }
}

UOutletBoundary::UOutletBoundary(const UOutletBoundary& ptf):fixedValueFvPatchField<vector>(ptf)
{
    if (false)
    {
        Info<<"construct ramU sha1: 2cb7a0c4345a61b88d329a26e832e625c53cdd92"
            " as copy\n";
    }
}

UOutletBoundary::UOutletBoundary(const UOutletBoundary& ptf,const DimensionedField<vector, volMesh>& iF):fixedValueFvPatchField<vector>(ptf, iF)
{
    if (false)
    {
        Info<<"construct ramU sha1: 2cb7a0c4345a61b88d329a26e832e625c53cdd92 "
            "as copy/DimensionedField\n";
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
UOutletBoundary::~UOutletBoundary()
{
    if (false)
    {
        Info<<"destroy ramU sha1: 2cb7a0c4345a61b88d329a26e832e625c53cdd92\n";
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void UOutletBoundary::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    if (false)
    {
        Info<<"updateCoeffs ramU sha1: 2cb7a0c4345a61b88d329a26e832e625c53cdd92\n";
    }
//begin code
    const scalar & time = this->db().time().value();   // Время
    const Foam::scalarIOList & scalarParameters = this->db().template lookupObject<scalarIOList>("scalarParameters");
    const Foam::vectorField & Cf = this->patch().Cf();
     // Получим текущий патч - границу
    const label ThisPatch = patch().index();
    vectorField & result = *this;
    // на граничной ячейки    
    tmp<vectorField> Field(this->patchInternalField());
    tmp<scalarField> rho = this->db().template lookupObject<volScalarField>("rho").boundaryField()[ThisPatch].patchInternalField();
    tmp<scalarField> p = this->db().template lookupObject<volScalarField>("p").boundaryField()[ThisPatch].patchInternalField();
    tmp<vectorField> U = this->db().template lookupObject<volVectorField>("U").boundaryField()[ThisPatch].patchInternalField();
    tmp<scalarField> gamma = this->db().template lookupObject<volScalarField>("gamma").boundaryField()[ThisPatch].patchInternalField();
    //-----------------------------------------------------
    const rhoReactionThermo & thermo = this->db().lookupObject<rhoReactionThermo>("thermophysicalProperties");
    rhoRDEThermo & Thermo = (rhoRDEThermo &)thermo;
    forAll(Cf, facei)
    {
        Foam::scalar MolWeight = Thermo.patchMolWeight(ThisPatch, facei);
        result[facei] = Uout(scalarParameters, time, MolWeight, p.ref()[facei], rho.ref()[facei], Field.ref()[facei], gamma.ref()[facei]);
    }
//this->operator==(vector(0,0,0));
// end code
    this->fixedValueFvPatchField<vector>::updateCoeffs();
}

} // End namespace Foam
