#include "TOutletBoundary.H"
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
// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //
//{{{ begin localCode
//}}} end localCode
// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //
extern "C"
{
    // unique function name that can be checked if the correct library version has been loaded
    void TOutletBoundary_1_0(bool load)
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
    TOutletBoundary
);
const char * const TOutletBoundary::SHA1sum = "66d4e0b38a18f2b7674332f6d021721a92f0e80e";

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
TOutletBoundary::TOutletBoundary(const fvPatch& p,const DimensionedField<scalar, volMesh>& iF):fixedValueFvPatchField<scalar>(p, iF)
{
    if (false)
    {
        Info<<"construct TOutletBoundary sha1: 66d4e0b38a18f2b7674332f6d021721a92f0e80e from patch/DimensionedField\n";
    }
}

TOutletBoundary::TOutletBoundary(const TOutletBoundary& ptf,const fvPatch& p,const DimensionedField<scalar, volMesh>& iF,const fvPatchFieldMapper& mapper):fixedValueFvPatchField<scalar>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct TOutletBoundary from patch/DimensionedField/mapper" << endl;
    }
}

TOutletBoundary::TOutletBoundary(const fvPatch& p,const DimensionedField<scalar, volMesh>& iF,const dictionary& dict):fixedValueFvPatchField<scalar>(p, iF, dict)
{
    if (false)
    {
        Info<<"construct TOutletBoundary from patch/dictionary" << endl;
    }
}

TOutletBoundary::TOutletBoundary(const TOutletBoundary& ptf):fixedValueFvPatchField<scalar>(ptf)
{
    if (false)
    {
        Info<<"construct TOutletBoundary as copy" << endl;;
    }
}

TOutletBoundary::TOutletBoundary(const TOutletBoundary& ptf, const DimensionedField<scalar, volMesh>& iF):
fixedValueFvPatchField<scalar>(ptf, iF)
{
    if (false)
    {
        Info<<"construct TOutletBoundary as copy/DimensionedField\n";
    }
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
TOutletBoundary::~TOutletBoundary()
{
    if (false)
    {
        Info << "destroy TOutletBoundary" << endl;
    }
}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void TOutletBoundary::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    if (false)
    {
        Info << "updateCoeffs TOutletBoundary" << endl;
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
    tmp<scalarField> T = this->db().template lookupObject<volScalarField>("T").boundaryField()[ThisPatch].patchInternalField();
    //-----------------------------------------------------
    const rhoReactionThermo & thermo = this->db().lookupObject<rhoReactionThermo>("thermophysicalProperties");
    rhoRDEThermo & Thermo = (rhoRDEThermo &)thermo;
    forAll(Cf, facei)
    {
        const Foam::scalar MolWeight = Thermo.patchMolWeight(ThisPatch, facei);
        result[facei] = ::Tout(scalarParameters, time, MolWeight, p.ref()[facei], rho.ref()[facei], U.ref()[facei], gamma.ref()[facei]);
    }
//this->operator==(300.0);
// end code
    this->fixedValueFvPatchField<scalar>::updateCoeffs();
}
} // End namespace Foam

/*
    Пример граничной ячейки:
    const label patchi = this->patch().template index();
    const scalarField TempIntF(Temp.boundaryField()[patchi].patchInternalField());
    const volScalarField Temp = db().lookupObject<volScalarField>(TempName_);
    или так:
    const label patchi = this->patch().template index();
    const volScalarField Temp = this->db().template lookupObject<volScalarField>(TempName_);
    Пример на границе:
    const scalarField & rho = this->patch().lookupPatchField<volScalarField, scalar>("rho");
*/