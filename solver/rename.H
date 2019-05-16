
Foam::volScalarField & PSI = const_cast<Foam::volScalarField &> (psi);
RenamePatchVolume(PSI, "cyclic", "calculated");
RenamePatchVolume(PSI, "processorCyclic", "calculated");

RenamePatchVolume(rho, "cyclic", "timeWallCyclic");
RenamePatchVolume(rho, "processorCyclic", "timeWallProcessorCyclicVolume");
rho.correctBoundaryConditions();

RenamePatchVolume(e, "cyclic", "timeWallCyclic");
RenamePatchVolume(e, "processorCyclic", "timeWallProcessorCyclicVolume");
e.correctBoundaryConditions();


RenamePatchVolume(p, "processorCyclic", "timeWallProcessorCyclicVolume");
p.correctBoundaryConditions();

RenamePatchVolume(T, "processorCyclic", "timeWallProcessorCyclicVolume");
T.correctBoundaryConditions();

RenamePatchVolume(MolWeight, "processorCyclic", "timeWallProcessorCyclicVolume");
MolWeight.correctBoundaryConditions();

RenamePatchVolume(Induction, "processorCyclic", "timeWallProcessorCyclicVolume");
Induction.correctBoundaryConditions();

RenamePatchVolume(U, "processorCyclic", "timeWallProcessorCyclicVolume");
U.correctBoundaryConditions();



/*
RenamePatchSurface(Rhof, "processorCyclic", "timeWallProcessorCyclicSurface");
RenamePatchSurface(Uf, "processorCyclic", "timeWallProcessorCyclicSurface");
RenamePatchSurface(Pf, "processorCyclic", "timeWallProcessorCyclicSurface");
RenamePatchSurface(Ef, "processorCyclic", "timeWallProcessorCyclicSurface");
RenamePatchSurface(Ethermodynamicalf, "processorCyclic", "timeWallProcessorCyclicSurface");
RenamePatchSurface(MolWeightf, "processorCyclic", "timeWallProcessorCyclicSurface");
RenamePatchSurface(Inductionf, "processorCyclic", "timeWallProcessorCyclicSurface");
*/

Info << e.boundaryField().types() << endl;
Info << Ef.boundaryField().types() << endl;
