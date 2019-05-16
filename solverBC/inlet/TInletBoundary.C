#include "TInletBoundary.H"
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
    void TInletBoundary_1_0(bool load)
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
    TInletBoundary
);
const char * const TInletBoundary::SHA1sum = "66d4e0b38a18f2b7674332f6d021721a92f0e80e";

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
TInletBoundary::TInletBoundary(const fvPatch& p,const DimensionedField<scalar, volMesh>& iF):fixedValueFvPatchField<scalar>(p, iF)
{
    if (false)
    {
        Info<<"construct TInletBoundary sha1: 66d4e0b38a18f2b7674332f6d021721a92f0e80e from patch/DimensionedField\n";
    }
}

TInletBoundary::TInletBoundary(const TInletBoundary& ptf,const fvPatch& p,const DimensionedField<scalar, volMesh>& iF,const fvPatchFieldMapper& mapper):fixedValueFvPatchField<scalar>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct TInletBoundary from patch/DimensionedField/mapper" << endl;
    }
}

TInletBoundary::TInletBoundary(const fvPatch& p,const DimensionedField<scalar, volMesh>& iF,const dictionary& dict):fixedValueFvPatchField<scalar>(p, iF, dict)
{
    if (false)
    {
        Info<<"construct TInletBoundary from patch/dictionary" << endl;
    }
}

TInletBoundary::TInletBoundary(const TInletBoundary& ptf):fixedValueFvPatchField<scalar>(ptf)
{
    if (false)
    {
        Info<<"construct TInletBoundary as copy" << endl;;
    }
}

TInletBoundary::TInletBoundary(const TInletBoundary& ptf, const DimensionedField<scalar, volMesh>& iF):fixedValueFvPatchField<scalar>(ptf, iF)
{
    if (false)
    {
        Info<<"construct TInletBoundary as copy/DimensionedField\n";
    }
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
TInletBoundary::~TInletBoundary()
{
    if (false)
    {
        Info << "destroy TInletBoundary" << endl;
    }
}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void TInletBoundary::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    if (false)
    {
        Info << "updateCoeffs TInletBoundary" << endl;
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
        result[facei] = Tin(scalarParameters, time, MolWeight, p.ref()[facei], rho.ref()[facei], U.ref()[facei], gamma.ref()[facei]);
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