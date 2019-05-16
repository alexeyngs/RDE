{
        fileName outputFile(runTime.path()/runTime.timeName());
        outputFile += "#";
        OFstream os(outputFile);
        os << "count = " << long(iter) << endl;
        os << "\t Induction  InductionZmin  InductionZmax Zmin \t Zmax \t T \t MolWeight" << endl;
        forAll(Induction, i)
        {
            os.width(3); os << i;
            os.width(12); os << Induction[i];
            os.width(12); os << InductionZmin[i];
            os.width(12); os << InductionZmax[i];
            os.width(12); os << Zmin[i];
            os.width(12); os << Zmax[i];
            os.width(12); os << T[i];
            os.width(12); os << MolWeight[i];
            os << endl;
        }
        os << endl;
        os << "==============================================================================================================" << endl;
        os << "\t Induction  Inductionf  InductionZmin  InductionZmax Zmin \t Zmax \t T \t MolWeight" << endl;
        forAll(Induction.boundaryField(), patchi) forAll(Induction.boundaryField()[patchi], i)
        {
            os.width(3); os << patchi; os.width(3); os << "-"; os.width(3); os << i;
            os.width(12); os << Induction.boundaryField()[patchi][i];
            os.width(12); os << Inductionf.boundaryField()[patchi][i];
            os.width(12); os << InductionZmin.boundaryField()[patchi][i];
            os.width(12); os << InductionZmax.boundaryField()[patchi][i];
            os.width(12); os << Zmin.boundaryField()[patchi][i];
            os.width(12); os << Zmax.boundaryField()[patchi][i];
            os.width(12); os << T.boundaryField()[patchi][i];
            os.width(12); os << MolWeight.boundaryField()[patchi][i];
            os << endl;
        }
        os << endl;
        forAll(ContactCoordinate, i) os << i << " * " << ContactCoordinate[i] << endl;
}
