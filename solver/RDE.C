#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "pointPatch.H"
#include "scalarIOList.H"

#include "thermo/rhoRDEThermo.H"
void Stitch(Foam::Time & runTime, Foam::polyMesh & mesh, const Foam::word masterPatchName, const Foam::word slavePatchName);
void ReadBreakDown(const Foam::dictionary & dict);
void BreakDown(Foam::rhoRDEThermo & Thermo, const Foam::volScalarField & Psi,
    const Foam::volScalarField & Rho, const Foam::volVectorField & U, const Foam::volScalarField & P, const Foam::volScalarField & E, const Foam::volScalarField & MolWeight, const Foam::volScalarField & Induction,
    const Foam::surfaceScalarField & Rholf, const Foam::surfaceVectorField & Ulf, const Foam::surfaceScalarField & Plf, const Foam::surfaceScalarField & Elf, const Foam::surfaceScalarField & MolWeightlf, const Foam::surfaceScalarField & Inductionlf, const Foam::surfaceScalarField & clf, const Foam::surfaceScalarField & gammalf,
    Foam::surfaceScalarField & Rhof, Foam::surfaceVectorField & Uf, Foam::surfaceScalarField & Pf, Foam::surfaceScalarField & Ef, Foam::surfaceScalarField & Ethermodynamicalf, Foam::surfaceScalarField & MolWeightf, Foam::surfaceScalarField & Inductionf);
inline Foam::scalar ABS(Foam::scalar x)
{
    return x < 0 ? -x : x;
}
const Foam::scalar Epsilon = Foam::SMALL*10.0;
void cellface(const Foam::fvMesh & mesh, Foam::label celli)
{
    
    const Foam::labelList & cFaces = mesh.cells()[celli];
    const Foam::point & centrevector = mesh.cellCentres()[celli];
    Foam::Info << "celli = " << celli << " - ";
    forAll(cFaces, cFacei)
    {
        Foam::label facei = cFaces[cFacei];
        const Foam::point & facevector = mesh.faceCentres()[facei];
        Foam::scalar delta = mag(facevector - centrevector);
        Info << facei << "|";
    }
    Foam::Info << endl;
}

template<typename TYPE>
inline void RenamePatchVolume(Foam::GeometricField<TYPE, Foam::fvPatchField, Foam::volMesh> & Field, word oldtype, word newtype)
{
    typename Foam::GeometricField<TYPE, Foam::fvPatchField, Foam::volMesh>::Boundary & Boundarys = Field.boundaryFieldRef();
    forAll(Boundarys, i)
    {
        if(Boundarys[i].type() == oldtype)
        {
            Foam::fvPatchField<TYPE> & PatchField = Boundarys[i];
            Boundarys.set(i, fvPatchField<TYPE>::New(newtype, oldtype, Field.mesh().boundary()[i], PatchField.internalField()));
        }
    }
};

template<typename TYPE>
inline void RenamePatchSurface(Foam::GeometricField<TYPE, Foam::fvsPatchField, Foam::surfaceMesh> & Field, word oldtype, word newtype)
{
    typename Foam::GeometricField<TYPE, Foam::fvsPatchField, Foam::surfaceMesh>::Boundary & Boundarys = Field.boundaryFieldRef();
    forAll(Boundarys, i)
    {
        if(Boundarys[i].type() == oldtype)
        {
            Foam::fvsPatchField<TYPE> & PatchField = Boundarys[i];
            Boundarys.set(i, fvsPatchField<TYPE>::New(newtype, oldtype, Field.mesh().boundary()[i], PatchField.internalField()));
        }
    }
};
Foam::volScalarField * FieldRho = nullptr;
Foam::volScalarField * FieldRhoE = nullptr;
bool TimeWall = true;
//----------------------------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #include "postProcess.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createBoundary.H"
    #include "createFields.H"
    #include "createFace.H"
    #include "rename.H"
    #include "createTimeControls.H"
    // Настройки с файла controlDict
    Foam::scalar CameraLength = readScalar(runTime.controlDict().lookup("CameraLength"));
    Foam::scalar CameraDiameter = readScalar(runTime.controlDict().lookup("CameraDiameter"));
    //dimensionedScalar p0("p0", dimensionSet(1, -1, -2, 0, 0, 0, 0), readScalar(runTime.controlDict().lookup("p0")));
    dimensionedScalar p0("p0", dimensionSet(1, -1, -2, 0, 0, 0, 0), 101325.0);
    bool UseChemistry = readBool(runTime.controlDict().lookup("UseChemistry"));
    //----------------------------------------------------------------------------------------------------------------
    #include "createContactFields.H"
    // Courant numbers used to adjust the time-step
    Courant = (mag(U) + c)*runTime.deltaT() / hdelta;
    Foam::scalar CoNum = max(Courant).value();
    // Инициализация переключателя на периодические условия
    Foam::scalar & timeOnePeriod = scalarParameters[0]; // Время, за которое детонационная волна пройдет один период
    timeOnePeriod = CameraDiameter * Foam::constant::mathematical::pi / thermo.Dcj();
    TimeWall = runTime.value() < timeOnePeriod * 0.9;
    ReadBreakDown(runTime.controlDict());
    //----------------------------------------------------------------------------------------------------------------
    Foam::Info << "My rank in MPI Communicator is " << Pstream::myProcNo() << " and master rank " << Pstream::masterNo() << Foam::endl;
    Foam::Info << "timeOnePeriod = " << timeOnePeriod << Foam::endl;
    Foam::label iter = runTime.startTimeIndex();
    Foam::Info << "Start iterator = " << iter << Foam::endl;
    Foam::Info << "Count cells = " << MolWeight.size() << Foam::endl;
    Foam::Info << "Size long double = " << sizeof(long double) << Foam::endl;
    Foam::Info << "CameraLength = " << CameraLength << Foam::endl;
    Foam::Info << "CameraDiameter = " << CameraDiameter << Foam::endl;
    Foam::Info << "UseChemistry = " << UseChemistry << Foam::endl;
    Foam::Info << "Starting time LOOP!" << Foam::endl;
    Foam::Info << "==========================================================================================================================" << endl << endl;
    while (runTime.run())
    {
        BreakDown(thermo, psi, rho, U, p, e, MolWeight, Induction, Rholf, Ulf, Plf, Elf, MolWeightlf, Inductionlf, clf, gammalf, Rhof, Uf, Pf, Ef, Ethermodynamicalf, MolWeightf, Inductionf);
        #include "updateFace.H"
        #include "readTimeControls.H"
        #include "setDeltaT.H"
        //<<<<<<<<<<<<<<<<<<<<<<<<
        Courant = (mag(U) + c)*runTime.deltaT() / hdelta;
        CoNum = max(Courant).value();
        Foam::Info << "Courant Number = " << CoNum << Foam::endl;
        Foam::Info << "Delta T: " << runTime.deltaTValue() << Foam::endl;
        //>>>>>>>>>>>>>>>>>>>>>>>>
        // Переключение от жесткой стенки к периодическим условиям
        TimeWall = runTime.value() < timeOnePeriod * 0.9;
        runTime++;
        iter++;
        Foam::Info << "Time = " << runTime.timeName() << nl << Foam::endl;
        // Уравнение массы
        solve
        (
            fvm::ddt(rho) + fvc::div(RhoUf)
        );

        // Уравнение импульса
        solve
        (
            fvm::ddt(rhoU) + fvc::div(RhoUUf) + fvc::div(PPf)
        );

        // Вычисление скорости
        U.ref() = rhoU() / rho();
        U.correctBoundaryConditions();
        rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();

        // Уравнение энергии
        solve
        (
            fvm::ddt(rhoE) + fvc::div(RhoUHf)
        );
        e = rhoE/rho - 0.5*magSqr(U);
        e.correctBoundaryConditions();
        rhoE.boundaryFieldRef() == rho.boundaryField() * (e.boundaryField() + 0.5*magSqr(U.boundaryField()));

        if(UseChemistry) thermo.chemistry(CameraLength);
        
        // Уравнение доли периода индукции Induction
        solve
        (
            fvm::ddt(rhoInduction) + fvc::div(RhoUInductionf) - FInduction
        );
        Induction.ref() = rhoInduction() / rho();
        forAll(Induction, i)
        {
            if(Induction[i] < 0.0) Induction[i] = 0.0;
            if(ABS(Zmax[i] - CameraLength) < Epsilon) Induction[i] = 1.0;
        }
        Induction.correctBoundaryConditions();
        rhoInduction.ref() = rho*Induction;
        rhoInduction.boundaryFieldRef() == rho.boundaryField()*Induction.boundaryField();

        // Уравнение молярной массы MolWeight
        solve
        (
            fvm::ddt(rhoMolWeight) + fvc::div(RhoUMolWeightf) - FMolWeight
        );
        // Подготовка молярной массы
        MolWeight.ref() = rhoMolWeight() / rho();
        thermo.CorrectMinMax();
        // Контактная граница
        // Интерполяци по термодинамической части энергии
        // Ethermodynamicalf -> Field
        UpdateContactField(mesh, Zmin, Zmax, Phimin, Phimax, InductionZmin, InductionZmax, FieldZmin, FieldZmax, PPhimin, PPhimax, Inductionf, Ethermodynamicalf, Pf);

        // Интерполяция по молярной массе
        // MolWeightf -> Field
//        UpdateContactField(mesh, Zmin, Zmax, Phimin, Phimax, InductionZmin, InductionZmax, FieldZmin, FieldZmax, PPhimin, PPhimax, Inductionf, MolWeightf, Pf);

        if(runTime.value() > timeOnePeriod/4.0) Contact(thermo, Induction, MolWeight, U, Zmin, Zmax, InductionZmin, InductionZmax, FieldZmin, FieldZmax, PPhimin, PPhimax, ContactCoordinateData, ContactCoordinateMesh, runTime.deltaTValue());
        MolWeight.correctBoundaryConditions();
        // Конец создания молярной массы
        rhoMolWeight.ref() = rho*MolWeight;
        rhoMolWeight.boundaryFieldRef() == rho.boundaryField()*MolWeight.boundaryField();
        //----------------------------------------------------------------------
        // Нахождение температуры
        thermo.correct();
        // Исправление ошибок по температуре
        thermo.CorrectErrors(UseChemistry, rho, Rhof, Pf);

        // Correct pressure
        p.ref() = rho() / psi();
        p.correctBoundaryConditions();
        rho.boundaryFieldRef() == psi.boundaryField() * p.boundaryField();

        g = -rho*U.component(vector::Z);
        impuls = p + rho*U.component(vector::Z)*U.component(vector::Z) - p0;
        runTime.write();
        Foam::Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl << Foam::endl;
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//#include "file/FieldFile.H"
//#include "file/ContactFile.H"
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        forAll(p, celli)
        {
            if(T[celli] < 0.0 || p[celli] < 0.0 || rho[celli] < 0.0)
            {
                Foam::scalar PP = p[celli];
                Foam::scalar TT = T[celli];
                Foam::scalar rhorho = rho[celli];
                Foam::Info << "CONTROL ERROR! celli = " << celli << " P = " << PP << " , T = " << TT << " , rho = " << rhorho << Foam::endl;
            }
        }
    }
    Foam::Info << "End" << Foam::endl;
    return 0;
}

template<typename TScalar> TScalar SQRT(TScalar X);
template<typename TScalar> TScalar POWER(TScalar X, TScalar Y); // X^Y
//====================================================================================================================
#include <cmath>
template<> Foam::scalar SQRT<Foam::scalar>(Foam::scalar X)
{
    return Foam::sqrt(X);
};
template<> Foam::scalar POWER<Foam::scalar>(Foam::scalar X, Foam::scalar Y) // X^Y
{
	return Foam::pow(X,Y);
};


/*
"includePath": [
                "${workspaceFolder}/**",
                "/opt/openfoam6/src/thermophysicalModels/reactionThermo/lnInclude/",
                "/opt/openfoam6/src/finiteVolume/cfdTools/general/include",
                "/opt/openfoam6/src/finiteVolume/finiteVolume/ddtSchemes/localEulerDdtScheme",
                "/opt/openfoam6/src/finiteVolume/lnInclude",
                "/opt/openfoam6/src/OpenFOAM/lnInclude",
                "/opt/openfoam6/src/TurbulenceModels/compressible/lnInclude",
                "/opt/openfoam6/src/TurbulenceModels/turbulenceModels/lnInclude",
                "/opt/openfoam6/src/thermophysicalModels/basic/lnInclude",
                "/opt/openfoam6/src/transportModels/compressible/lnInclude",
                "/opt/openfoam6/src/dynamicMesh/lnInclude",
                "/opt/openfoam6/src/mesh/blockMesh/lnInclude",
                "/opt/openfoam6/src/OSspecific/POSIX/lnInclude",
                "/opt/openfoam6/src/thermophysicalModels/reactionThermo/lnInclude",
                "/opt/openfoam6/src/thermophysicalModels/specie/lnInclude",
                "/usr/include/c++/7",
                "/usr/include/x86_64-linux-gnu/c++/7",
                "/opt/openfoam5/src/OSspecific/POSIX"
            ],
*/
