Info<< "Reading thermophysical properties\n" << endl;


autoPtr<Foam::rhoRDEThermo> pThermo
(
    rhoRDEThermo::New(mesh)
);

Foam::rhoRDEThermo & thermo = pThermo();

thermo.correct();

Foam::volScalarField & e = thermo.he();
Foam::volScalarField & p = thermo.p();
Foam::volScalarField & T = thermo.T();
const Foam::volScalarField & psi = thermo.psi();
const Foam::volScalarField & mu = thermo.mu();
Foam::volScalarField & MolWeight = thermo.MolWeight();
Foam::volScalarField & FMolWeight = thermo.FMolWeight();
Foam::volScalarField & Induction = thermo.Induction();
Foam::volScalarField & FInduction = thermo.FInduction();

autoPtr<compressible::turbulenceModel> turbulence;
autoPtr<volScalarField> muEffPtr;
autoPtr<volTensorField> tauMCPtr;
autoPtr<surfaceScalarField> sigmaDotUPtr;


Foam::volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Foam::volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

Foam::volVectorField rhoU
(
    IOobject
    (
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U
);

Foam::volScalarField rhoE
(
    IOobject
    (
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*(e + 0.5*magSqr(U))
);

Foam::volScalarField rhoMolWeight
(
    IOobject
    (
        "rhoMolWeight",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*MolWeight
);

Foam::volScalarField rhoInduction
(
    IOobject
    (
        "rhoInduction",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*Induction
);


Foam::volScalarField rPsi
(
    "rPsi",
    1.0/psi
);

Foam::volScalarField gamma
(
    "gamma",
    thermo.Cp()/thermo.Cv()
);

Foam::volScalarField c
(
    "c",
    sqrt(gamma * rPsi)
);

// OpenFOAM delta coeffs can be negative for coupled patches
Foam::volScalarField hdelta
(
    IOobject
    (
        "hdelta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    1.0
);

//<<<<<<<<<<<<<<<<<<<<<<<<<<<
    forAll(mesh.cells(), celli) // По всем ячейкам
    {
        Foam::scalar lenght = GREAT; // Максимально возможный размер
        const Foam::labelList & cFaces = mesh.cells()[celli];
        const Foam::point & centrevector = mesh.cellCentres()[celli];
        forAll(cFaces, cFacei)
        {
            Foam::label facei = cFaces[cFacei];
            const Foam::point & facevector = mesh.faceCentres()[facei];
            Foam::scalar delta = mag(facevector - centrevector);
            if(lenght > delta) lenght = delta;
        }
        hdelta[celli] = lenght;
    }
//>>>>>>>>>>>>>>>>>>>>>>>>>>>

Foam::volScalarField Courant
(
    "Courant",
    (mag(U) + c)*runTime.deltaT() / hdelta
);
//==========================================================================================================

Foam::volScalarField g
(
    IOobject
    (
        "g",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U.component(vector::Z)
);

Foam::volScalarField impuls
(
    IOobject
    (
        "impuls",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    p + rho*U.component(vector::Z)*U.component(vector::Z)
);