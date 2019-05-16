template<typename TScalar, typename TYPE>
inline TYPE gradient(const TScalar dH, const TYPE ff1, const TYPE ff2, const TYPE ff3)
// dH - расстояние от центра до грани - половина размера ячейки
{
    TYPE outresult = (ff3 - ff1)/(dH*2.0);
    if(mag(outresult) > mag(ff3 - ff2)/dH) outresult = (ff3 - ff2)/dH;
    if(mag(outresult) > mag(ff2 - ff1)/dH) outresult = (ff2 - ff1)/dH;
    return outresult;
};

template<typename TScalar, typename TYPE>
inline TYPE gradient(const TScalar L1, const TScalar L2, const TYPE ff1, const TYPE ff2, const TYPE ff3)
// L1 - размер левой ячейки
// L2 - размер правой ячейки
{
    TYPE outresult = (ff3 - ff1)/((L1 + L2) / 2.0);
    if(mag(outresult) > mag(ff3 - ff2) / (L2 / 2.0)) outresult = (ff3 - ff2) / (L2 / 2.0);
    if(mag(outresult) > mag(ff2 - ff1) / (L1 / 2.0)) outresult = (ff2 - ff1) / (L1 / 2.0);
    return outresult;
};

inline Foam::label OpposideFaceLabel(const fvMesh & mesh, Foam::label celli, Foam::label face)
{
    // Координата текущей грани
    const point & facevector = mesh.faceCentres()[celli];
    // Список всех граней
    const Foam::labelList & cFaces = mesh.cells()[celli];
    Foam::label result = 0;
    Foam::scalar maxdelta = 0.0;
    forAll(cFaces, cFacei)
    {
        Foam::label facei = cFaces[cFacei]; // получаем индекс сетки
        const Foam::point & p = mesh.faceCentres()[facei];  // центр ячейки с определеннаи индексом
        Foam::scalar delta = mag(p - facevector);
        if(delta > maxdelta)
        {
            delta = maxdelta;
            result = facei;
        }
    }
    return result;
};

template<typename TYPE>
inline TYPE Face(const Foam::GeometricField<TYPE, Foam::fvsPatchField, Foam::surfaceMesh> & Field, Foam::label index)
{
    if(index < Field.size()) return Field[index];
    else
    {
        index -= Field.size();
        forAll(Field.boundaryField(), patchi)
        {
            if(index >= Field.boundaryField()[patchi].size())
            {
                index -= Field.boundaryField()[patchi].size();
            }else
            {
                return Field.boundaryField()[patchi][index];
            }
        }
    }
    return Field[0];
};
//=================================================================================================================
template<typename TScalar, typename TVector, typename SFType>
inline void Kolgan
(
    rhoRDEThermo & Thermo,
    const GeometricField<TScalar, fvPatchField, volMesh> & Psi,
    const SFType & Sf,
    const GeometricField<TScalar, fvPatchField, volMesh> & Rho,
    const GeometricField<TVector, fvPatchField, volMesh> & U,
    const GeometricField<TScalar, fvPatchField, volMesh> & p,
    const GeometricField<TScalar, fvPatchField, volMesh> & E,
    const GeometricField<TScalar, fvPatchField, volMesh> & MolWeight,
    const GeometricField<TScalar, fvPatchField, volMesh> & Induction,
    const GeometricField<TScalar, fvPatchField, volMesh> & c,
    const GeometricField<TScalar, fvPatchField, volMesh> & gamma,
    const Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & Rholf,
    const Foam::GeometricField<TVector, Foam::fvsPatchField, Foam::surfaceMesh> & Ulf,
    const Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & Plf,
    const Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & Elf,
    const Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & MolWeightlf,
    const Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & Inductionlf,
    const Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & clf,
    const Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & gammalf,
    Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & Rhof,
    Foam::GeometricField<TVector, Foam::fvsPatchField, Foam::surfaceMesh> & Uf,
    Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & Pf,
    Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & Ef,
    Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & Ethermodynamicalf,
    Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & MolWeightf,
    Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & Inductionf
)
{
    const fvMesh & mesh = Rho.mesh();
    const labelUList & P = mesh.owner();
    const labelUList & N = mesh.neighbour();
    //-------------------------------------------------------------------------------
    const tmp<surfaceScalarField> & tlambdas = Rho.mesh().surfaceInterpolation::weights();
    if (surfaceInterpolation::debug)
    {
        InfoInFunction << "Interpolating " << Rho.type() << " " << Rho.name() << " from cells to faces without explicit correction" << endl;
    }

    const surfaceScalarField & lambdas = tlambdas();
    const scalarField & lambda = lambdas;

    const Field<TScalar> & Rho_ = Rho;
    const Field<TVector> & U_ = U;
    const Field<TScalar> & p_ = p;
    const Field<TScalar> & E_ = E;
    const Field<TScalar> & MolWeight_ = MolWeight;
    const Field<TScalar> & Induction_ = Induction;
    const Field<TScalar> & gamma_ = gamma;
    const Field<TScalar> & c_ = c;
    
    Field<TScalar> & Rhof_ = Rhof.primitiveFieldRef();
    Field<TVector> & Uf_ = Uf.primitiveFieldRef();
    Field<TScalar> & Pf_ = Pf.primitiveFieldRef();
    Field<TScalar> & Ef_ = Ef.primitiveFieldRef();
    Field<TScalar> & Ethermodynamicalf_ = Ethermodynamicalf.primitiveFieldRef();
    Field<TScalar> & MolWeightf_ = MolWeightf.primitiveFieldRef();
    Field<TScalar> & Inductionf_ = Inductionf.primitiveFieldRef();

    const GeometricField<TVector, fvsPatchField, surfaceMesh> & _Sf = mesh.Sf();
    const GeometricField<TScalar, fvsPatchField, surfaceMesh> & _magSf = mesh.magSf();
    const GeometricField<TVector, fvPatchField, volMesh> & CellCentres = mesh.C();
    const GeometricField<TVector, fvsPatchField, surfaceMesh> & FaceCentres = mesh.Cf();

    for (label fi = 0; fi < P.size(); fi++)
    {
        TVector Normal = _Sf[fi] / _magSf[fi];
		// Ячейка №1 - P
        TScalar dP = mag(FaceCentres[fi] - CellCentres[P[fi]]);
        TScalar rho1 = Rho_[P[fi]];
		TScalar p1 = p_[P[fi]];
        TVector u1 = U_[P[fi]];
        TScalar gamma1 = gamma_[P[fi]];
        TScalar E1 = E_[P[fi]];
        TScalar MolWeight1 = MolWeight_[P[fi]];
        TScalar Induction1 = Induction_[P[fi]];
        // Поправка Колгана для ячейки №1 - P
        label OpposideP = OpposideFaceLabel(mesh, P[fi], fi); // Индекс противоположной грани
        rho1 += gradient(dP, Face(Rholf, OpposideP), rho1, Rholf[fi])*dP;
        p1 += gradient(dP, Face(Plf, OpposideP), p1, Plf[fi])*dP;
        u1 += gradient(dP, Face(Ulf, OpposideP), u1, Ulf[fi])*dP;
        E1 += gradient(dP, Face(Elf, OpposideP), E1, Elf[fi])*dP;
        MolWeight1 += gradient(dP, Face(MolWeightlf, OpposideP), MolWeight1, MolWeightlf[fi])*dP;
        Induction1 += gradient(dP, Face(Inductionlf, OpposideP), Induction1, Inductionlf[fi])*dP;
        E1 -= Thermo.Ethermodynamical(p1, rho1, MolWeight1);
        // Разложение скорости
        TScalar U1 = u1 & Normal;		// Проекция скорости на нормаль - нормальная компонента скорости в ячейке 1
		TVector T1 = u1 - (Normal * U1);	// Тангенсальная компонента скорости в ячейке 1
        TScalar c1 = c_[P[fi]];
		// Ячейка №2 - N
        TScalar dN =  mag(CellCentres[N[fi]] - FaceCentres[fi]);
        TScalar rho2 = Rho_[N[fi]];
		TScalar p2 = p_[N[fi]];
        TVector u2 = U_[N[fi]];
        TScalar gamma2 = gamma_[N[fi]];
        TScalar E2 = E_[N[fi]];
        TScalar MolWeight2 = MolWeight_[N[fi]];
        TScalar Induction2 = Induction_[N[fi]];
        // Поправка Колгана для ячейки №2 - N
        label OpposideN = OpposideFaceLabel(mesh, N[fi], fi); // Индекс противоположной грани
        rho2 -= gradient(dN, Rholf[fi], rho2, Face(Rholf, OpposideN))*dN;
        p2 -= gradient(dN, Plf[fi], p2, Face(Plf, OpposideN))*dN;
        u2 -= gradient(dN, Ulf[fi], u2, Face(Ulf, OpposideN))*dN;
        E2 -= gradient(dN, Elf[fi], E2, Face(Elf, OpposideN))*dN;
        MolWeight2 -= gradient(dN, MolWeightlf[fi], MolWeight2, Face(MolWeightlf, OpposideN))*dN;
        Induction2 -= gradient(dN, Inductionlf[fi], Induction2, Face(Inductionlf, OpposideN))*dN;
        E2 -= Thermo.Ethermodynamical(p2, rho2, MolWeight2);
		// Разложение скорости
		TScalar U2 = u2 & Normal;		// Проекция скорости на нормаль - нормальная компонента скорости в ячейке 2
		TVector T2 = u2 - (Normal * U2);	// Тангенсальная компонента скорости в ячейке 2
        TScalar c2 = c_[N[fi]];
        //----------------------------------------------------
		// out
		TScalar PP;
		TScalar RR;
		TScalar UU;
		TScalar Ustar;                      // скорость контактного разрыва
		TScalar W = 0.0;
		int sln = SLN<TScalar>(p1, rho1, U1, gamma1, c1, p2, rho2, U2, gamma2, c2, PP, RR, UU, Ustar, W);
        if(sln == -5 || sln == +5) Info << "Vacuum" << endl;
		Rhof_[fi] = RR;
        Pf_[fi] = PP;
		if(sln < 0)
		{
            Uf_[fi] = T1 + (Normal * UU);
        	MolWeightf_[fi] = MolWeight1;
        	Inductionf_[fi] = Induction1;
            Ethermodynamicalf_[fi] = Thermo.Ethermodynamical(PP, RR, MolWeight1);
			Ef_[fi] = E1 + Ethermodynamicalf_[fi];
		}else if(sln > 0)
		{
            Uf_[fi] = T2 + (Normal * UU);
        	MolWeightf_[fi] = MolWeight2;
        	Inductionf_[fi] = Induction2;
            Ethermodynamicalf_[fi] = Thermo.Ethermodynamical(PP, RR, MolWeight2);
			Ef_[fi] = E2 + Ethermodynamicalf_[fi];
		}else
        {
            Info << "ERORR!" << endl;
            break;	//Error
        }
    }
    // Interpolate across coupled patches using given lambdas
    typename GeometricField<TScalar, fvsPatchField, surfaceMesh>::Boundary & Rhobf = Rhof.boundaryFieldRef();
    typename GeometricField<TVector, fvsPatchField, surfaceMesh>::Boundary & Ubf = Uf.boundaryFieldRef();
    typename GeometricField<TScalar, fvsPatchField, surfaceMesh>::Boundary & Pbf = Pf.boundaryFieldRef();
    typename GeometricField<TScalar, fvsPatchField, surfaceMesh>::Boundary & Ebf = Ef.boundaryFieldRef();
    typename GeometricField<TScalar, fvsPatchField, surfaceMesh>::Boundary & Ethermodynamicalbf = Ethermodynamicalf.boundaryFieldRef();
    typename GeometricField<TScalar, fvsPatchField, surfaceMesh>::Boundary & MolWeightbf = MolWeightf.boundaryFieldRef();
    typename GeometricField<TScalar, fvsPatchField, surfaceMesh>::Boundary & Inductionbf = Inductionf.boundaryFieldRef();

    forAll(mesh.boundaryMesh(), ipatch)
    {
        word BCtype = mesh.boundaryMesh().types()[ipatch];
        const fvsPatchScalarField & pLambda = lambdas.boundaryField()[ipatch];
        const typename SFType::Patch & pSf = Sf.boundaryField()[ipatch];

        fvsPatchField<TScalar> & Rhosf = Rhobf[ipatch];
        fvsPatchField<TVector> & Usf = Ubf[ipatch];
        fvsPatchField<TScalar> & Psf = Pbf[ipatch];
        fvsPatchField<TScalar> & Esf = Ebf[ipatch];
        fvsPatchField<TScalar> & Ethermodynamicalsf = Ethermodynamicalbf[ipatch];
        fvsPatchField<TScalar> & MolWeightsf = MolWeightbf[ipatch];
        fvsPatchField<TScalar> & Inductionsf = Inductionbf[ipatch];

        // 2 в 1:
        // 1 - межпроцессорное взаимодействие processor
        // 2 - периодическое граничное условие processorPatch
        bool UseProcessorBounary = ((TimeWall && BCtype == "processor") || (!TimeWall && isA<processorFvPatch>(mesh.boundary()[ipatch])));
        // граничное условия для двух областей, но обе области нахидятся под одним процессором
        bool UsePatchBounary = (BCtype != "processor") && (BCtype != "cyclic") && (BCtype != "processorCyclic");
        // Если межпроцессорное взаимодействие
        if (Pstream::parRun() && mesh.boundaryMesh()[ipatch].coupled() && UseProcessorBounary)
        {
            const fvsPatchField<TVector> & Sf_ = mesh.Sf().boundaryField()[ipatch];
            const fvsPatchField<TScalar> & magSf_ = mesh.magSf().boundaryField()[ipatch];
            const fvPatchField<TVector> & CellCentres_ = mesh.C().boundaryField()[ipatch];
            const fvsPatchField<TVector> & FaceCentres_ = mesh.Cf().boundaryField()[ipatch];
            const labelList & boundaryFaces = mesh.boundaryMesh()[ipatch].meshEdges();  // Список индексов граней на границе            

            // Для конкретного патча P - Internal - Полоса ячеек этого процессора
            Field<TScalar> RhoInternal = Rho.boundaryField()[ipatch].patchInternalField();
            Field<TVector> UInternal = U.boundaryField()[ipatch].patchInternalField();
            Field<TScalar> PInternal = p.boundaryField()[ipatch].patchInternalField();
            Field<TScalar> EInternal = E.boundaryField()[ipatch].patchInternalField();
            Field<TScalar> MolWeightInternal = MolWeight.boundaryField()[ipatch].patchInternalField();
            Field<TScalar> InductionInternal = Induction.boundaryField()[ipatch].patchInternalField();
            Field<TScalar> gammaInternal = gamma.boundaryField()[ipatch].patchInternalField();
            Field<TScalar> cInternal = c.boundaryField()[ipatch].patchInternalField();
            Field<TScalar> deltaInternal = mag(FaceCentres_ - CellCentres_);
            const labelList & InternalCells = mesh.C().boundaryField()[ipatch].patch().patch().faceCells(); // Список индексов ячеек на границе(пограничные ячейки)
            // Полоски(патчи) грани от линейной интерполяции(lf) - slf
            const fvsPatchField<TScalar> & Rhoslf = Rholf.boundaryField()[ipatch];
            const fvsPatchField<TVector> & Uslf = Ulf.boundaryField()[ipatch];
            const fvsPatchField<TScalar> & Pslf = Plf.boundaryField()[ipatch];
            const fvsPatchField<TScalar> & Eslf = Elf.boundaryField()[ipatch];
            const fvsPatchField<TScalar> & MolWeightslf = MolWeightlf.boundaryField()[ipatch];
            const fvsPatchField<TScalar> & Inductionslf = Inductionlf.boundaryField()[ipatch];
            const fvsPatchField<TScalar> & gammaslf = gammalf.boundaryField()[ipatch];
            const fvsPatchField<TScalar> & cslf = clf.boundaryField()[ipatch];
            // Заполнение поправок
            forAll(Rho.boundaryField()[ipatch], i)
            {
                TScalar delta = deltaInternal[i];
                label Edge = boundaryFaces[i]; // Индекс грани
                label Opposide = OpposideFaceLabel(mesh, InternalCells[i], Edge); // Индекс противоположной грани
                RhoInternal[i] += gradient(delta, Face(Rholf, Opposide), RhoInternal[i], Rhoslf[i])*delta;
                UInternal[i] += gradient(delta, Face(Ulf, Opposide), UInternal[i], Uslf[i])*delta;
                PInternal[i] += gradient(delta, Face(Plf, Opposide), PInternal[i], Pslf[i])*delta;
                EInternal[i] += gradient(delta, Face(Elf, Opposide), EInternal[i], Eslf[i])*delta;
                MolWeightInternal[i] += gradient(delta, Face(MolWeightlf, Opposide), MolWeightInternal[i], MolWeightslf[i])*delta;
                InductionInternal[i] += gradient(delta, Face(Inductionlf, Opposide), InductionInternal[i], Inductionslf[i])*delta;
                gammaInternal[i] += gradient(delta, Face(gammalf, Opposide), gammaInternal[i], gammaslf[i])*delta;
                cInternal[i] += gradient(delta, Face(clf, Opposide), cInternal[i], cslf[i])*delta;
            }
            // Для конкретного патча N - Neighbour - Полоса ячеек соседнего процессора
            Field<TScalar> RhoNeighbour;
            Field<TVector> UNeighbour;
            Field<TScalar> PNeighbour;
            Field<TScalar> ENeighbour;
            Field<TScalar> MolWeightNeighbour;
            Field<TScalar> InductionNeighbour;
            Field<TScalar> gammaNeighbour;
            Field<TScalar> cNeighbour;
            Field<TScalar> deltaNeighbour;
            // Получение идентификатора процессора(себя и соседнего)
            const Foam::processorPolyPatch & ProcessorPolyPatch = refCast<const Foam::processorPolyPatch>(mesh.boundaryMesh()[ipatch]);
            // Обмен данных с соседним процессором(отдаем свои - Internal и получаем от соседа - Neighbour)
            MPIexchange<TScalar, TVector>
            (
                RhoInternal, UInternal, PInternal, EInternal, MolWeightInternal, InductionInternal, gammaInternal, cInternal, deltaInternal,
                RhoNeighbour, UNeighbour, PNeighbour, ENeighbour, MolWeightNeighbour, InductionNeighbour, gammaNeighbour, cNeighbour, deltaNeighbour,
                ProcessorPolyPatch
            );
            forAll(Rho.boundaryField()[ipatch], i)
            {
                TVector Normal = Sf_[i] / magSf_[i];
                // Ячейка №1 - P
                TScalar rho1 = RhoInternal[i];
                TScalar p1 = PInternal[i];
                TVector u1 = UInternal[i];
                TScalar gamma1 = gammaInternal[i];
                TScalar E1 = EInternal[i];
                TScalar MolWeight1 = MolWeightInternal[i];
                TScalar Induction1 = InductionInternal[i];
                E1 -= Thermo.Ethermodynamical(p1, rho1, MolWeight1);
                // Разложение скорости
                TScalar U1 = u1 & Normal;		// Проекция скорости на нормаль - нормальная компонента скорости в ячейке 1
                TVector T1 = u1 - (Normal * U1);	// Тангенсальная компонента скорости в ячейке 1
                TScalar c1 = cInternal[i];
                // Ячейка №2 - N
                TScalar rho2 = RhoNeighbour[i];
                TScalar p2 = PNeighbour[i];
                TVector u2 = UNeighbour[i];
                TScalar gamma2 = gammaNeighbour[i];
                TScalar E2 = ENeighbour[i];
                TScalar MolWeight2 = MolWeightNeighbour[i];
                TScalar Induction2 = InductionNeighbour[i];
                // Поправка Колгана для ячейки №2 - N
                E2 -= Thermo.Ethermodynamical(p2, rho2, MolWeight2);
                // Разложение скорости
                TScalar U2 = u2 & Normal;		// Проекция скорости на нормаль - нормальная компонента скорости в ячейке 2
                TVector T2 = u2 - (Normal * U2);	// Тангенсальная компонента скорости в ячейке 2
                TScalar c2 = cNeighbour[i];
                //----------------------------------------------------
                // out
                TScalar PP;
                TScalar RR;
                TScalar UU;
                TScalar Ustar;                      // скорость контактного разрыва
                TScalar W = 0.0;
                int sln = SLN<TScalar>(p1, rho1, U1, gamma1, c1, p2, rho2, U2, gamma2, c2, PP, RR, UU, Ustar, W);
                if(sln == -5 || sln == +5) Info << "Vacuum" << endl;
                Rhosf[i] = RR;
                Psf[i] = PP;
                if(sln < 0)
                {
                    Usf[i] = T1 + (Normal * UU);
                    MolWeightsf[i] = MolWeight1;
                    Inductionsf[i] = Induction1;
                    Ethermodynamicalsf[i] = Thermo.Ethermodynamical(PP, RR, MolWeight1);
                    Esf[i] = E1 + Ethermodynamicalsf[i];
                }else if(sln > 0)
                {
                    Usf[i] = T2 + (Normal * UU);
                    MolWeightsf[i] = MolWeight2;
                    Inductionsf[i] = Induction2;
                    Ethermodynamicalsf[i] = Thermo.Ethermodynamical(PP, RR, MolWeight2);
                    Esf[i] = E2 + Ethermodynamicalsf[i];
                }else
                {
                    Info << "ERORR!" << endl;
                    break;	//Error
                }
            }
            RhoInternal.clear();
            UInternal.clear();
            PInternal.clear();
            EInternal.clear();
            MolWeightInternal.clear();
            InductionInternal.clear();
            gammaInternal.clear();
            cInternal.clear();
            deltaInternal.clear();
            //------------------------------------------------------------
            RhoNeighbour.clear();
            UNeighbour.clear();
            PNeighbour.clear();
            ENeighbour.clear();
            MolWeightNeighbour.clear();
            InductionNeighbour.clear();
            gammaNeighbour.clear();
            cNeighbour.clear();
            deltaNeighbour.clear();
        }else if (mesh.boundaryMesh()[ipatch].coupled() && UsePatchBounary) // PatchIndex
        {
            goto patch;
        }
        else
        {
            patch:
            forAll(Rho.boundaryField()[ipatch], i)
            {
                Rhosf[i] = Rho.boundaryField()[ipatch][i];
                Usf[i] = U.boundaryField()[ipatch][i];
                Psf[i] = p.boundaryField()[ipatch][i];
                Esf[i] = E.boundaryField()[ipatch][i];
                MolWeightsf[i] = MolWeight.boundaryField()[ipatch][i];
                Inductionsf[i] = Induction.boundaryField()[ipatch][i];
                Ethermodynamicalsf[i] = Thermo.Ethermodynamical(p.boundaryField()[ipatch][i], Rho.boundaryField()[ipatch][i], MolWeight.boundaryField()[ipatch][i]);
            }
        }
    }
    tlambdas.clear();
}
