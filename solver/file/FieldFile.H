{
        fileName outputFile(runTime.path()/runTime.timeName());
        OFstream os(outputFile);

        os << "count cells: " << p.size() << "  Step number: " <<  long(iter) << "   Courant: " << CoNum << endl;
        forAll(ContactCoordinate, i) 
        os << i << " - " << ContactCoordinate[i] << " - " << ContactTEST[i] << endl;
        os << "======================================================================================" << endl;
        os << "\t\tP \t rho \t MolWeight \t   T \t\te \t\t rhoE \t\t gamma \t\t FMolWeight \t\t U" << endl;
        forAll(p, i)
        {
            os.width(6); os << i;
            os.width(12); os << p[i]/101325.0;
            os.width(12); os << rho[i];
            os.width(12); os << MolWeight[i];
            os.width(12); os << T[i];
            os.width(16); os << e[i];
            os.width(16); os << rhoE[i];
            os.width(16); os << gamma[i];
            os.width(12); os << FMolWeight[i];
            os.width(12); os << U[i];
            os << endl;
        }
        os << "--------------------------------------------------------------------------------------------------------------" << endl;
        os << "\t\tPf \trhof  \tMolWeightf \t\t Ef \t\t Hf \t\t RhoUf \t\t RhoUHf\tRhoUMolWeightf \t\t Ethermodynamicalf \t\t Uf" << endl;
        forAll(Pf, i)
        {
            os.width(6); os << i;
            os.width(12); os << Pf[i]/101325.0;
            os.width(12); os << Rhof[i];
            os.width(12); os << MolWeightf[i];
            os.width(16); os << Ef[i];
            os.width(16); os << Hf[i];
            os.width(16); os << RhoUf[i];
            os.width(16); os << RhoUHf[i];
            os.width(14); os << RhoUMolWeightf[i];
            os.width(14); os << Ethermodynamicalf[i];
            os.width(14); os << Uf[i];
            os << endl;
        }
        /*
        os << "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000" << endl;
        os << "\t\tPf \trhof  \tMolWeightf \t Ef \t\t Hf \t\t RhoUHf \t Uf" << endl;
        forAll(p.boundaryField(), patchi) forAll(p.boundaryField()[patchi], i)
        {
            os.width(3); os << patchi; os.width(3); os << "-"; os.width(3); os << i;
            os.width(9); os << p.boundaryField()[patchi][i]/101325.0;
            os.width(12); os << rho.boundaryField()[patchi][i];
            os.width(12); os << MolWeight.boundaryField()[patchi][i];
            os.width(16); os << e.boundaryField()[patchi][i];
            os.width(16); os << Hf.boundaryField()[patchi][i];
            os.width(12); os << RhoUHf.boundaryField()[patchi][i];
            os.width(12); os << U.boundaryField()[patchi][i];
            os << endl;
        }*/
        os << endl;
        os << "==============================================================================================================" << endl;
        // Ребра
        const cellList & cells = mesh.cells();
        scalarField hmax(cells.size());
/*        forAll(mesh.cells(),celli)
        {
            const labelList & cFaces = mesh.cells()[celli];
            const point & centrevector = mesh.cellCentres()[celli];
            forAll(cFaces, cFacei)
            {
                label facei = cFaces[cFacei];
                const point & facevector = mesh.faceCentres()[facei];
                scalar delta = mag(facevector - centrevector);
                os << facei << "|";
            }
            os << endl;
        }
*/
}//getchar();