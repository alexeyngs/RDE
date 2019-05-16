template<typename TScalar, typename TVector, typename SFType>
inline void Accustic
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
    const tmp<surfaceScalarField> & tlambdas = mesh.surfaceInterpolation::weights();
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
		int sln = SLN_accustic<TScalar>(p1, rho1, U1, gamma1, c1, p2, rho2, U2, gamma2, c2, PP, RR, UU, Ustar, W);
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

//        if (Pstream::parRun() && mesh.boundaryMesh()[ipatch].coupled() && BCtype == "processor")//- if running parallel
        if (Pstream::parRun() && mesh.boundaryMesh()[ipatch].coupled() && isA<processorFvPatch>(mesh.boundary()[ipatch])) //- if running parallel
        {
            const fvsPatchField<TVector> & Sf_ = mesh.Sf().boundaryField()[ipatch];
            const fvsPatchField<TScalar> & magSf_ = mesh.magSf().boundaryField()[ipatch];
            const fvPatchField<TVector> & CellCentres_ = mesh.C().boundaryField()[ipatch];
            const fvsPatchField<TVector> & FaceCentres_ = mesh.Cf().boundaryField()[ipatch];

            // Для конкретного патча P - Internal - Полоса ячеек этого процессора
            tmp<Field<TScalar>> RhoInternal = Rho.boundaryField()[ipatch].patchInternalField();
            tmp<Field<TVector>> UInternal = U.boundaryField()[ipatch].patchInternalField();
            tmp<Field<TScalar>> PInternal = p.boundaryField()[ipatch].patchInternalField();
            tmp<Field<TScalar>> EInternal = E.boundaryField()[ipatch].patchInternalField();
            tmp<Field<TScalar>> MolWeightInternal = MolWeight.boundaryField()[ipatch].patchInternalField();
            tmp<Field<TScalar>> InductionInternal = Induction.boundaryField()[ipatch].patchInternalField();
            tmp<Field<TScalar>> gammaInternal = gamma.boundaryField()[ipatch].patchInternalField();
            tmp<Field<TScalar>> cInternal = c.boundaryField()[ipatch].patchInternalField();
            tmp<Field<TScalar>> deltaInternal = mag(FaceCentres_ - CellCentres_);

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
                int sln = SLN_accustic<TScalar>(p1, rho1, U1, gamma1, c1, p2, rho2, U2, gamma2, c2, PP, RR, UU, Ustar, W);
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
        }else if (mesh.boundaryMesh()[ipatch].coupled() && BCtype != "processor") // PatchIndex
        {
            const fvsPatchField<TVector> & Sf_ = mesh.Sf().boundaryField()[ipatch];
            const fvsPatchField<TScalar> & magSf_ = mesh.magSf().boundaryField()[ipatch];
            const fvPatchField<TVector> & CellCentres_ = mesh.C().boundaryField()[ipatch];
            const fvsPatchField<TVector> & FaceCentres_ = mesh.Cf().boundaryField()[ipatch];

            // Для конкретного патча поля для P
            tmp<Field<TScalar>> RhoInternal = Rho.boundaryField()[ipatch].patchInternalField();
            tmp<Field<TVector>> UInternal = U.boundaryField()[ipatch].patchInternalField();
            tmp<Field<TScalar>> PInternal = p.boundaryField()[ipatch].patchInternalField();
            tmp<Field<TScalar>> EInternal = E.boundaryField()[ipatch].patchInternalField();
            tmp<Field<TScalar>> MolWeightInternal = MolWeight.boundaryField()[ipatch].patchInternalField();
            tmp<Field<TScalar>> InductionInternal = Induction.boundaryField()[ipatch].patchInternalField();
            tmp<Field<TScalar>> gammaInternal = gamma.boundaryField()[ipatch].patchInternalField();
            tmp<Field<TScalar>> cInternal = c.boundaryField()[ipatch].patchInternalField();

            // Для конкретного патча параметры N
            tmp<Field<TScalar>> RhoNeighbour = Rho.boundaryField()[ipatch].patchNeighbourField();
            tmp<Field<TVector>> UNeighbour = U.boundaryField()[ipatch].patchNeighbourField();
            tmp<Field<TScalar>> PNeighbour = p.boundaryField()[ipatch].patchNeighbourField();            
            tmp<Field<TScalar>> ENeighbour = E.boundaryField()[ipatch].patchNeighbourField();            
            tmp<Field<TScalar>> MolWeightNeighbour = MolWeight.boundaryField()[ipatch].patchNeighbourField();            
            tmp<Field<TScalar>> InductionNeighbour = Induction.boundaryField()[ipatch].patchNeighbourField();            
            tmp<Field<TScalar>> gammaNeighbour = gamma.boundaryField()[ipatch].patchNeighbourField();
            tmp<Field<TScalar>> cNeighbour = c.boundaryField()[ipatch].patchNeighbourField();
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
                int sln = SLN_accustic<TScalar>(p1, rho1, U1, gamma1, c1, p2, rho2, U2, gamma2, c2, PP, RR, UU, Ustar, W);
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
