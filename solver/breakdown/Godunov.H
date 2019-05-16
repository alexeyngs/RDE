template<typename TScalar, typename TVector, typename SFType>
inline void Godunov
(
    Foam::rhoRDEThermo & Thermo,
    const Foam::GeometricField<TScalar, fvPatchField, volMesh> & Psi,
    const SFType & Sf,
    const Foam::GeometricField<TScalar, Foam::fvPatchField, Foam::volMesh> & Rho,
    const Foam::GeometricField<TVector, Foam::fvPatchField, Foam::volMesh> & U,
    const Foam::GeometricField<TScalar, Foam::fvPatchField, Foam::volMesh> & p,
    const Foam::GeometricField<TScalar, Foam::fvPatchField, Foam::volMesh> & E,
    const Foam::GeometricField<TScalar, Foam::fvPatchField, Foam::volMesh> & MolWeight,
    const Foam::GeometricField<TScalar, Foam::fvPatchField, Foam::volMesh> & Induction,
    const Foam::GeometricField<TScalar, Foam::fvPatchField, Foam::volMesh> & c,
    const Foam::GeometricField<TScalar, Foam::fvPatchField, Foam::volMesh> & gamma,
    Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & Rhof,
    Foam::GeometricField<TVector, Foam::fvsPatchField, Foam::surfaceMesh> & Uf,
    Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & Pf,
    Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & Ef,
    Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & Ethermodynamicalf,
    Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & MolWeightf,
    Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh> & Inductionf
)
{
    const Foam::fvMesh & mesh = Rho.mesh();
    const Foam::labelUList & P = mesh.owner();
    const Foam::labelUList & N = mesh.neighbour();
    //-------------------------------------------------------------------------------
    const Foam::tmp<surfaceScalarField> & tlambdas = mesh.surfaceInterpolation::weights();
    if (surfaceInterpolation::debug)
    {
        InfoInFunction << "Interpolating " << Rho.type() << " " << Rho.name() << " from cells to faces without explicit correction" << endl;
    }

    const Foam::surfaceScalarField & lambdas = tlambdas();
    const Foam::scalarField & lambda = lambdas;
    
    const Foam::Field<TScalar> & Rho_ = Rho;
    const Foam::Field<TVector> & U_ = U;
    const Foam::Field<TScalar> & p_ = p;
    const Foam::Field<TScalar> & E_ = E;
    const Foam::Field<TScalar> & MolWeight_ = MolWeight;
    const Foam::Field<TScalar> & Induction_ = Induction;
    const Foam::Field<TScalar> & gamma_ = gamma;
    const Foam::Field<TScalar> & c_ = c;
    
    Foam::Field<TScalar> & Rhof_ = Rhof.primitiveFieldRef();
    Foam::Field<TVector> & Uf_ = Uf.primitiveFieldRef();
    Foam::Field<TScalar> & Pf_ = Pf.primitiveFieldRef();
    Foam::Field<TScalar> & Ef_ = Ef.primitiveFieldRef();
    Foam::Field<TScalar> & Ethermodynamicalf_ = Ethermodynamicalf.primitiveFieldRef();
    Foam::Field<TScalar> & MolWeightf_ = MolWeightf.primitiveFieldRef();
    Foam::Field<TScalar> & Inductionf_ = Inductionf.primitiveFieldRef();
    const Foam::GeometricField<TVector, fvsPatchField, surfaceMesh> & _Sf = mesh.Sf();
    const Foam::GeometricField<TScalar, fvsPatchField, surfaceMesh> & _magSf = mesh.magSf();
    const Foam::GeometricField<TVector, fvPatchField, volMesh> & CellCentres = mesh.C();
    const Foam::GeometricField<TVector, fvsPatchField, surfaceMesh> & FaceCentres = mesh.Cf();

    for (Foam::label fi = 0; fi < P.size(); fi++)
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

    // Interpolate across coupled patches using given lamb  das
    typename Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh>::Boundary & Rhobf = Rhof.boundaryFieldRef();
    typename Foam::GeometricField<TVector, Foam::fvsPatchField, Foam::surfaceMesh>::Boundary & Ubf = Uf.boundaryFieldRef();
    typename Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh>::Boundary & Pbf = Pf.boundaryFieldRef();
    typename Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh>::Boundary & Ebf = Ef.boundaryFieldRef();
    typename Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh>::Boundary & Ethermodynamicalbf = Ethermodynamicalf.boundaryFieldRef();
    typename Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh>::Boundary & MolWeightbf = MolWeightf.boundaryFieldRef();
    typename Foam::GeometricField<TScalar, Foam::fvsPatchField, Foam::surfaceMesh>::Boundary & Inductionbf = Inductionf.boundaryFieldRef();

    forAll(mesh.boundaryMesh(), ipatch)
    {
        Foam::word BCtype = mesh.boundaryMesh().types()[ipatch];
        const Foam::fvsPatchScalarField & pLambda = lambdas.boundaryField()[ipatch];
        const typename SFType::Patch & pSf = Sf.boundaryField()[ipatch];

        Foam::fvsPatchField<TScalar> & Rhosf = Rhobf[ipatch];
        Foam::fvsPatchField<TVector> & Usf = Ubf[ipatch];
        Foam::fvsPatchField<TScalar> & Psf = Pbf[ipatch];
        Foam::fvsPatchField<TScalar> & Esf = Ebf[ipatch];
        Foam::fvsPatchField<TScalar> & Ethermodynamicalsf = Ethermodynamicalbf[ipatch];
        Foam::fvsPatchField<TScalar> & MolWeightsf = MolWeightbf[ipatch];
        Foam::fvsPatchField<TScalar> & Inductionsf = Inductionbf[ipatch];

        // 2 в 1:
        // 1 - межпроцессорное взаимодействие processor
        // 2 - периодическое граничное условие processorPatch
        bool UseProcessorBounary = ((TimeWall && BCtype == "processor") || (!TimeWall && isA<processorFvPatch>(mesh.boundary()[ipatch])));
        // граничное условия для двух областей, но обе области нахидятся под одним процессором
        bool UsePatchBounary = (BCtype != "processor") && (BCtype != "cyclic") && (BCtype != "processorCyclic");
        // Если межпроцессорное взаимодействие
        if (Pstream::parRun() && mesh.boundaryMesh()[ipatch].coupled() && UseProcessorBounary)
        {
            const Foam::fvsPatchField<TVector> & Sf_ = mesh.Sf().boundaryField()[ipatch];
            const Foam::fvsPatchField<TScalar> & magSf_ = mesh.magSf().boundaryField()[ipatch];
            const Foam::fvPatchField<TVector> & CellCentres_ = mesh.C().boundaryField()[ipatch];
            const Foam::fvsPatchField<TVector> & FaceCentres_ = mesh.Cf().boundaryField()[ipatch];

            // Для конкретного патча P - Internal - Полоса ячеек этого процессора
            Foam::tmp<Foam::Field<TScalar>> RhoInternal = Rho.boundaryField()[ipatch].patchInternalField();
            Foam::tmp<Foam::Field<TVector>> UInternal = U.boundaryField()[ipatch].patchInternalField();
            Foam::tmp<Foam::Field<TScalar>> PInternal = p.boundaryField()[ipatch].patchInternalField();
            Foam::tmp<Foam::Field<TScalar>> EInternal = E.boundaryField()[ipatch].patchInternalField();
            Foam::tmp<Foam::Field<TScalar>> MolWeightInternal = MolWeight.boundaryField()[ipatch].patchInternalField();
            Foam::tmp<Foam::Field<TScalar>> InductionInternal = Induction.boundaryField()[ipatch].patchInternalField();
            Foam::tmp<Foam::Field<TScalar>> gammaInternal = gamma.boundaryField()[ipatch].patchInternalField();
            Foam::tmp<Foam::Field<TScalar>> cInternal = c.boundaryField()[ipatch].patchInternalField();
            Foam::tmp<Foam::Field<TScalar>> deltaInternal = mag(FaceCentres_ - CellCentres_);

            // Для конкретного патча N - Neighbour - Полоса ячеек соседнего процессора
            Foam::Field<TScalar> RhoNeighbour;
            Foam::Field<TVector> UNeighbour;
            Foam::Field<TScalar> PNeighbour;
            Foam::Field<TScalar> ENeighbour;
            Foam::Field<TScalar> MolWeightNeighbour;
            Foam::Field<TScalar> InductionNeighbour;
            Foam::Field<TScalar> gammaNeighbour;
            Foam::Field<TScalar> cNeighbour;
            Foam::Field<TScalar> deltaNeighbour;
            
            const Foam::processorPolyPatch & ProcessorPolyPatch = refCast<const Foam::processorPolyPatch>(mesh.boundaryMesh()[ipatch]);
            MPIexchange<TScalar, TVector>
            (
                RhoInternal.ref(), UInternal.ref(), PInternal.ref(), EInternal.ref(), MolWeightInternal.ref(), InductionInternal.ref(), gammaInternal.ref(), cInternal.ref(), deltaInternal.ref(),
                RhoNeighbour, UNeighbour, PNeighbour, ENeighbour, MolWeightNeighbour, InductionNeighbour, gammaNeighbour, cNeighbour, deltaNeighbour,
                ProcessorPolyPatch
            );
            forAll(Rho.boundaryField()[ipatch], i)
            {
                TVector Normal = Sf_[i] / magSf_[i];
                // Ячейка №1 - P
                TScalar dP = deltaInternal.ref()[i];
                TScalar rho1 = RhoInternal.ref()[i];
                TScalar p1 = PInternal.ref()[i];
                TVector u1 = UInternal.ref()[i];
                TScalar gamma1 = gammaInternal.ref()[i];
                TScalar E1 = EInternal.ref()[i];
                TScalar MolWeight1 = MolWeightInternal.ref()[i];
                TScalar Induction1 = InductionInternal.ref()[i];
                E1 -= Thermo.Ethermodynamical(p1, rho1, MolWeight1);
                // Разложение скорости
                TScalar U1 = u1 & Normal;		// Проекция скорости на нормаль - нормальная компонента скорости в ячейке 1
                TVector T1 = u1 - (Normal * U1);	// Тангенсальная компонента скорости в ячейке 1
                TScalar c1 = cInternal.ref()[i];
                // Ячейка №2 - N
                TScalar dN = deltaNeighbour[i];
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
            const Foam::fvsPatchField<TVector> & Sf_ = mesh.Sf().boundaryField()[ipatch];
            const Foam::fvsPatchField<TScalar> & magSf_ = mesh.magSf().boundaryField()[ipatch];
            const Foam::fvPatchField<TVector> & CellCentres_ = mesh.C().boundaryField()[ipatch];
            const Foam::fvsPatchField<TVector> & FaceCentres_ = mesh.Cf().boundaryField()[ipatch];

            // Для конкретного патча поля для P
            Foam::tmp<Foam::Field<TScalar>> RhoInternal = Rho.boundaryField()[ipatch].patchInternalField();
            Foam::tmp<Foam::Field<TVector>> UInternal = U.boundaryField()[ipatch].patchInternalField();
            Foam::tmp<Foam::Field<TScalar>> PInternal = p.boundaryField()[ipatch].patchInternalField();
            Foam::tmp<Foam::Field<TScalar>> EInternal = E.boundaryField()[ipatch].patchInternalField();
            Foam::tmp<Foam::Field<TScalar>> MolWeightInternal = MolWeight.boundaryField()[ipatch].patchInternalField();
            Foam::tmp<Foam::Field<TScalar>> InductionInternal = Induction.boundaryField()[ipatch].patchInternalField();
            Foam::tmp<Foam::Field<TScalar>> gammaInternal = gamma.boundaryField()[ipatch].patchInternalField();
            Foam::tmp<Foam::Field<TScalar>> cInternal = c.boundaryField()[ipatch].patchInternalField();

            // Для конкретного патча параметры N
            Foam::tmp<Foam::Field<TScalar>> RhoNeighbour = Rho.boundaryField()[ipatch].patchNeighbourField();
            Foam::tmp<Foam::Field<TVector>> UNeighbour = U.boundaryField()[ipatch].patchNeighbourField();
            Foam::tmp<Foam::Field<TScalar>> PNeighbour = p.boundaryField()[ipatch].patchNeighbourField();            
            Foam::tmp<Foam::Field<TScalar>> ENeighbour = E.boundaryField()[ipatch].patchNeighbourField();            
            Foam::tmp<Foam::Field<TScalar>> MolWeightNeighbour = MolWeight.boundaryField()[ipatch].patchNeighbourField();            
            Foam::tmp<Foam::Field<TScalar>> InductionNeighbour = Induction.boundaryField()[ipatch].patchNeighbourField();            
            Foam::tmp<Foam::Field<TScalar>> gammaNeighbour = gamma.boundaryField()[ipatch].patchNeighbourField();
            Foam::tmp<Foam::Field<TScalar>> cNeighbour = c.boundaryField()[ipatch].patchNeighbourField();
            forAll(Rho.boundaryField()[ipatch], i)
            {
                TVector Normal = Sf_[i] / magSf_[i];
                // Ячейка №1 - P
                TScalar dP = mag(FaceCentres_[i] - CellCentres_[i]);
                TScalar rho1 = RhoInternal.ref()[i];
                TScalar p1 = PInternal.ref()[i];
                TVector u1 = UInternal.ref()[i];
                TScalar gamma1 = gammaInternal.ref()[i];
                TScalar E1 = EInternal.ref()[i];
                TScalar MolWeight1 = MolWeightInternal.ref()[i];
                TScalar Induction1 = InductionInternal.ref()[i];
                E1 -= Thermo.Ethermodynamical(p1, rho1, MolWeight1);
                // Разложение скорости
                TScalar U1 = u1 & Normal;		// Проекция скорости на нормаль - нормальная компонента скорости в ячейке 1
                TVector T1 = u1 - (Normal * U1);	// Тангенсальная компонента скорости в ячейке 1
                TScalar c1 = cInternal.ref()[i];
                // Ячейка №2 - N
                TScalar dN = mag(CellCentres_[i] - FaceCentres_[i]);
                TScalar rho2 = RhoNeighbour.ref()[i];
                TScalar p2 = PNeighbour.ref()[i];
                TVector u2 = UNeighbour.ref()[i];
                TScalar gamma2 = gammaNeighbour.ref()[i];
                TScalar E2 = ENeighbour.ref()[i];
                TScalar MolWeight2 = MolWeightNeighbour.ref()[i];
                TScalar Induction2 = InductionNeighbour.ref()[i];
                // Поправка Колгана для ячейки №2 - N
                E2 -= Thermo.Ethermodynamical(p2, rho2, MolWeight2);
                // Разложение скорости
                TScalar U2 = u2 & Normal;		// Проекция скорости на нормаль - нормальная компонента скорости в ячейке 2
                TVector T2 = u2 - (Normal * U2);	// Тангенсальная компонента скорости в ячейке 2
                TScalar c2 = cNeighbour.ref()[i];
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
        }
        else
        {
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
