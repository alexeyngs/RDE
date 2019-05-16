
Foam::volScalarField Zmin
(
    IOobject
    (
        "Zmin",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    CameraLength
);

Foam::volScalarField Zmax
(
    IOobject
    (
        "Zmax",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE 
    ),
    mesh,
    0.0
);

Foam::volScalarField Phimin
(
    IOobject
    (
        "Phimin",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    2.0*Foam::constant::mathematical::pi
);

Foam::volScalarField Phimax
(
    IOobject
    (
        "Phimax",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    0.0
);
//-------------------------------------------------------------------------------------

Foam::volScalarField InductionZmin
(
    IOobject
    (
        "InductionZmin",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    0.0
);

Foam::volScalarField InductionZmax
(
    IOobject
    (
        "InductionZmax",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    0.0
);

Foam::volScalarField FieldZmin
(
    IOobject
    (
        "FieldZmin",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    0.0
);

Foam::volScalarField FieldZmax
(
    IOobject
    (
        "FieldZmax",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    0.0
);

Foam::volScalarField PPhimin
(
    IOobject
    (
        "PPhimin",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    0.0
);

Foam::volScalarField PPhimax
(
    IOobject
    (
        "PPhimax",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    0.0
);

void GetContactEdgeSize(const Foam::fvMesh & mesh, const Foam::scalar CameraLength, Foam::vectorField & ContactCoordinateMesh, Foam::scalarField & ContactCoordinate);
Foam::vectorField ContactCoordinateMesh;    // Сетка с координатами
Foam::scalarField ContactCoordinateData;    // Сетка с данными
// Получение координат и данных
GetContactEdgeSize(mesh, CameraLength, ContactCoordinateMesh, ContactCoordinateData);
//---------------------------------------------------------------------------------------------------------------------------
void Contact(Foam::rhoRDEThermo & thermo, const Foam::volScalarField & Induction, Foam::volScalarField & MolWeight, const Foam::volVectorField & U, const volScalarField & Zmin, const Foam::volScalarField & Zmax,
const Foam::volScalarField & InductionZmin, const Foam::volScalarField & InductionZmax, const Foam::volScalarField & FieldZmin, const Foam::volScalarField & FieldZmax, const Foam::volScalarField & FieldPhimin, const Foam::volScalarField & FieldPhimax,
Foam::scalarField & ContactCoordinateField, const Foam::vectorField & ContactCoordinateMesh, const Foam::scalar dt);

void UpdateContactMesh(Foam::volScalarField & Zmin, Foam::volScalarField & Zmax, Foam::volScalarField & Phimin, Foam::volScalarField & Phimax);

void UpdateContactField(const Foam::fvMesh & mesh, const Foam::volScalarField & Zmin, const Foam::volScalarField & Zmax, const Foam::volScalarField & Phimin, const Foam::volScalarField & Phimax,
Foam::volScalarField & FieldZmin1, Foam::volScalarField & FieldZmax1, Foam::volScalarField & FieldZmin2, Foam::volScalarField & FieldZmax2, 
Foam::volScalarField & FieldPhimin, Foam::volScalarField & FieldPhimax, Foam::surfaceScalarField & FieldZf1, Foam::surfaceScalarField & FieldZf2, Foam::surfaceScalarField & FieldPhif);

UpdateContactMesh(Zmin, Zmax, Phimin, Phimax);
