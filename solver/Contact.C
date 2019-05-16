#include "IOmanip.H"
#include "fvCFD.H"
#include "pTraits.H"
#include "Field.H"
#include "boundaryMesh.H"
#include "thermo/rhoRDEThermo.H"

inline const Foam::label GetLabelBounary(const Foam::fvMesh & mesh, const Foam::vectorField & ContactCoordinateMesh, const Foam::label celli)
{
    // Точка (X,Y) для конкретной ячейки
    const Foam::point Position(mesh.cellCentres()[celli].x(), mesh.cellCentres()[celli].y(), 0.0);
    Foam::label result = -1;
    Foam::scalar lenght = Foam::GREAT;
    forAll(ContactCoordinateMesh, i) // По всем ячейкам
    {
        Foam::point P(ContactCoordinateMesh[i].x(), ContactCoordinateMesh[i].y(), 0.0);
        Foam::scalar delta = Foam::mag(P - Position);
        if(lenght > delta)
        {
            lenght = delta;
            result = i;
        }
    }
    return result;
}

inline Foam::scalar ABS(Foam::scalar x)
{
    return x < 0 ? -x : x;
}

void UpdateContactMesh(Foam::volScalarField & Zmin, Foam::volScalarField & Zmax, Foam::volScalarField & Phimin, Foam::volScalarField & Phimax)
{
    const Foam::fvMesh & mesh = Zmin.mesh();
    const Foam::labelUList & owner = mesh.owner();
    const Foam::labelUList & neighbour = mesh.neighbour();
    forAll(owner, facei)
    {
        //---------------------------------------------------------
        const Foam::point & Position = mesh.faceCentres()[facei];
        scalar Z = Position.z();
        scalar phi = Foam::acos(Position.x() / Foam::sqrt(Position.x() * Position.x() + Position.y() * Position.y()));
        //---------------------------------------------------------
        if(Zmin[owner[facei]] > Z) Zmin[owner[facei]] = Z;
        if(Zmax[owner[facei]] < Z) Zmax[owner[facei]] = Z;
        if(Phimin[owner[facei]] > phi) Phimin[owner[facei]] = phi;
        if(Phimax[owner[facei]] < phi) Phimax[owner[facei]] = phi;
        //---------------------------------------------------------
        if(Zmin[neighbour[facei]] > Z) Zmin[neighbour[facei]] = Z;
        if(Zmax[neighbour[facei]] < Z) Zmax[neighbour[facei]] = Z;
        if(Phimin[neighbour[facei]] > phi) Phimin[neighbour[facei]] = phi;
        if(Phimax[neighbour[facei]] < phi) Phimax[neighbour[facei]] = phi;
        //---------------------------------------------------------
    }
    forAll(mesh.boundaryMesh(), patchi)
    {
        const labelUList & pFaceCells = mesh.boundary()[patchi].faceCells();
        const vectorField & faceCentres = mesh.Cf().boundaryField()[patchi];
        if(faceCentres.empty()) continue;
        forAll(mesh.boundaryMesh()[patchi], facei)
        {
            //---------------------------------------------------------
            const Foam::vector & Position = faceCentres[facei];
            scalar Z = Position.z();
            scalar phi = Foam::acos(Position.x() / Foam::sqrt(Position.x() * Position.x() + Position.y() * Position.y()));
            //---------------------------------------------------------
            if(Zmin[pFaceCells[facei]] > Z) Zmin[pFaceCells[facei]] = Z;
            if(Zmax[pFaceCells[facei]] < Z) Zmax[pFaceCells[facei]] = Z;
            if(Phimin[pFaceCells[facei]] > phi) Phimin[pFaceCells[facei]] = phi;
            if(Phimax[pFaceCells[facei]] < phi) Phimax[pFaceCells[facei]] = phi;
        }
    }
}

void UpdateContactField(const Foam::fvMesh & mesh, const Foam::volScalarField & Zmin, const Foam::volScalarField & Zmax, const Foam::volScalarField & Phimin, const Foam::volScalarField & Phimax,
Foam::volScalarField & FieldZmin1, Foam::volScalarField & FieldZmax1, Foam::volScalarField & FieldZmin2, Foam::volScalarField & FieldZmax2, 
Foam::volScalarField & FieldPhimin, Foam::volScalarField & FieldPhimax, Foam::surfaceScalarField & FieldZf1, Foam::surfaceScalarField & FieldZf2, Foam::surfaceScalarField & FieldPhif)
{
    const Foam::scalar Epsilon = Foam::SMALL*10.0;
    const Foam::labelUList & owner = mesh.owner();
    const Foam::labelUList & neighbour = mesh.neighbour();
    forAll(owner, facei)
    {
        //---------------------------------------------------------
        const Foam::point & Position = mesh.faceCentres()[facei];
        Foam::scalar Z = Position.z();
        Foam::scalar phi = Foam::acos(Position.x() / Foam::sqrt(Position.x() * Position.x() + Position.y() * Position.y()));
        //---------------------------------------------------------
        if(ABS(Zmin[owner[facei]] - Z) < Epsilon) FieldZmin1[owner[facei]] = FieldZf1[facei];
        if(ABS(Zmax[owner[facei]] - Z) < Epsilon) FieldZmax1[owner[facei]] = FieldZf1[facei];
        if(ABS(Zmin[owner[facei]] - Z) < Epsilon) FieldZmin2[owner[facei]] = FieldZf2[facei];
        if(ABS(Zmax[owner[facei]] - Z) < Epsilon) FieldZmax2[owner[facei]] = FieldZf2[facei];
        if(ABS(Phimin[owner[facei]] - phi) < Epsilon) FieldPhimin[owner[facei]] = FieldPhif[facei];
        if(ABS(Phimax[owner[facei]] - phi) < Epsilon) FieldPhimax[owner[facei]] = FieldPhif[facei];
        //---------------------------------------------------------
        if(ABS(Zmin[neighbour[facei]] - Z) < Epsilon) FieldZmin1[neighbour[facei]] = FieldZf1[facei];
        if(ABS(Zmax[neighbour[facei]] - Z) < Epsilon) FieldZmax1[neighbour[facei]] = FieldZf1[facei];
        if(ABS(Zmin[neighbour[facei]] - Z) < Epsilon) FieldZmin2[neighbour[facei]] = FieldZf2[facei];
        if(ABS(Zmax[neighbour[facei]] - Z) < Epsilon) FieldZmax2[neighbour[facei]] = FieldZf2[facei];
        if(ABS(Phimin[neighbour[facei]] - phi) < Epsilon) FieldPhimin[neighbour[facei]] = FieldPhif[facei];
        if(ABS(Phimax[neighbour[facei]] - phi) < Epsilon) FieldPhimax[neighbour[facei]] = FieldPhif[facei];
        //---------------------------------------------------------
    }
    forAll(mesh.boundaryMesh(), patchi)
    {
        const Foam::labelUList & pFaceCells = mesh.boundary()[patchi].faceCells();
        const Foam::vectorField & faceCentres = mesh.Cf().boundaryField()[patchi];
        if(faceCentres.empty()) continue;
        const Foam::fvsPatchField<Foam::scalar> & BoundFieldZf1 = FieldZf1.boundaryField()[patchi];
        const Foam::fvsPatchField<Foam::scalar> & BoundFieldZf2 = FieldZf2.boundaryField()[patchi];
        const Foam::fvsPatchField<Foam::scalar> & BoundFieldPhif = FieldPhif.boundaryField()[patchi];
        forAll(mesh.boundaryMesh()[patchi], facei)
        {
            //---------------------------------------------------------
            const Foam::vector & Position = faceCentres[facei];
            Foam::scalar Z = Position.z();
            Foam::scalar phi = Foam::acos(Position.x() / Foam::sqrt(Position.x() * Position.x() + Position.y() * Position.y()));
            //---------------------------------------------------------
            if(ABS(Zmin[pFaceCells[facei]] - Z) < Epsilon) FieldZmin1[pFaceCells[facei]] = BoundFieldZf1[facei];
            if(ABS(Zmax[pFaceCells[facei]] - Z) < Epsilon) FieldZmax1[pFaceCells[facei]] = BoundFieldZf1[facei];
            if(ABS(Zmin[pFaceCells[facei]] - Z) < Epsilon) FieldZmin2[pFaceCells[facei]] = BoundFieldZf2[facei];
            if(ABS(Zmax[pFaceCells[facei]] - Z) < Epsilon) FieldZmax2[pFaceCells[facei]] = BoundFieldZf2[facei];
            if(ABS(Phimin[pFaceCells[facei]] - phi) < Epsilon) FieldPhimin[pFaceCells[facei]] = BoundFieldPhif[facei];
            if(ABS(Phimax[pFaceCells[facei]] - phi) < Epsilon) FieldPhimax[pFaceCells[facei]] = BoundFieldPhif[facei];
        }
    }
}

inline bool IsContact(const Foam::scalar FieldZmin, const Foam::scalar FieldZmax, const Foam::scalar FieldPhimin, const Foam::scalar FieldPhimax, const Foam::vector U)
{
//    Foam::scalar Uphi = Foam::acos(U.x() / Foam::sqrt(U.x() * U.x() + U.y() * U.y()));
//    Foam::scalar Uphi = U.x();
//    bool result = (FieldZmin < Foam::SMALL) && (FieldZmax > Foam::SMALL) && (Uphi < 120.0);    // Критерий через скорость
//    bool result = (FieldZmin < Foam::SMALL) && (FieldZmax > Foam::SMALL) && (FieldPhimax/FieldPhimin < 1.2); // Критерий через давление
    bool result = (FieldZmin < Foam::SMALL) && (FieldZmax > Foam::SMALL);
    return result;
}

void Contact(Foam::rhoRDEThermo & thermo, const Foam::volScalarField & Induction, Foam::volScalarField & MolWeight, const Foam::volVectorField & U, const volScalarField & Zmin, const Foam::volScalarField & Zmax,
const Foam::volScalarField & InductionZmin, const Foam::volScalarField & InductionZmax, const Foam::volScalarField & FieldZmin, const Foam::volScalarField & FieldZmax, const Foam::volScalarField & FieldPhimin, const Foam::volScalarField & FieldPhimax,
Foam::scalarField & ContactCoordinateField, const Foam::vectorField & ContactCoordinateMesh, const Foam::scalar dt)
{
//if(ContactCoordinateField.empty()) return;
    const Foam::fvMesh & mesh = MolWeight.mesh();
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    forAll(MolWeight, i) if(IsContact(InductionZmin[i], InductionZmax[i], FieldPhimin[i], FieldPhimax[i], U[i]))
    {
        // Получение индекса поля плоского "блинчика"(ContactCoordinateField), в котором хранятся значения контактных координат
        // index-i сетки камеры(mesh) => index-j сетки "блинчика"(ContactCoordinateMesh)
        const Foam::label j = GetLabelBounary(mesh, ContactCoordinateMesh, i);
        // Получение ссылки на значение контактной координаты для i ячейки mesh и j ячейки ContactCoordinateField
        Foam::scalar & ContactCoordinate = ContactCoordinateField[j];
        //--------------------------------------------------------------------
        if(ContactCoordinate < Zmin[i]) ContactCoordinate = Zmin[i];
        if(ContactCoordinate > Zmax[i]) ContactCoordinate = Zmax[i];
        ContactCoordinate += dt*U[i].z();
        // alpha - доля остаточной горячей смеси, уменьшается до 0
        Foam::scalar alpha = (ContactCoordinate - Zmin[i]) / (Zmax[i] - Zmin[i]);
        if(alpha < 0.0) alpha = 0.0;
        if(alpha > 1.0) alpha = 1.0;
        // Интерполяция
        Foam::scalar Ethermodynamical = thermo.Ethermodynamical(i)*alpha + FieldZmax[i]*(1.0-alpha);
        MolWeight[i] = thermo.GetmolWeight(i, Ethermodynamical);
//none        MolWeight[i] = MolWeight[i]*alpha + FieldZmax[i]*(1.0-alpha);
    }
}


#include <vector>
#include <utility>
// first - Z
// second - Mesh
typedef std::pair<Foam::scalar, Foam::vectorField> CellCandidate;
//-----------------------------------------------------------------------------------------------
// Функция для заполнения сетки(ContactCoordinateMesh) и данные(ContactCoordinateField) для контактной границы
// Всю сетку режет на блины(или полоски) при Z=const,
// а затем выбирают блин с максимальным размером
void GetContactEdgeSize(const Foam::fvMesh & mesh, const Foam::scalar CameraLength, Foam::vectorField & ContactCoordinateMesh, Foam::scalarField & ContactCoordinateField)
{

    std::vector<CellCandidate> DATA;
    // резка на блины
    forAll(mesh.Cf(), i)
    {
        Foam::vector V = mesh.Cf()[i];
        bool added = false;
        for(CellCandidate & elem : DATA) if(ABS(V.z() - elem.first) < Foam::SMALL)
        {
            // Увеличиваем размер блина при определенным значения Z
            elem.second.append(V);
            added = true;
            break;
        }
        // Создаем новый блин
        if(!added)
        {
            CellCandidate cell;
            cell.first = V.z();
            cell.second.append(V);
            DATA.push_back(cell);
        }
    }
    // Выбираем блин с максимальным размером
    size_t CurrentZ = 0;
    for(size_t i = 0; i < DATA.size(); i++)
    if(DATA[CurrentZ].second.size() < DATA[i].second.size())
    {
        CurrentZ = i;
    }
    // Сетка с координатами
    ContactCoordinateMesh = DATA[CurrentZ].second;
    // Сетка с данными
    ContactCoordinateField.resize(ContactCoordinateMesh.size());
    forAll(ContactCoordinateField, i) ContactCoordinateField[i] = CameraLength;
    DATA.clear();
    return;

    //---------------------------------------------
    Foam::label Patch = mesh.boundaryMesh().findPatchID("inlet");
    if(Patch < 0) return;
    // Сетка с данными
    const size_t size = mesh.boundaryMesh()[Patch].size();
    ContactCoordinateField.resize(size);
    forAll(ContactCoordinateField, i) ContactCoordinateField[i] = CameraLength;

    // Сетка с координатами
    ContactCoordinateMesh = mesh.Cf().boundaryField()[Patch];
}