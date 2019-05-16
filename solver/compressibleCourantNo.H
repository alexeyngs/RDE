/*
Global
    compressibleCourantNo
Description
    Calculates the mean and maximum wave speed based Courant Numbers.
*/
Foam::scalar CoNum = 0.0;
Foam::scalar meanCoNum = 0.0;
if (mesh.nInternalFaces())
{
    surfaceScalarField amaxSfbyDelta
    (
        mesh.surfaceInterpolation::deltaCoeffs()*amaxSf
    );
    CoNum = max(amaxSfbyDelta/mesh.magSf()).value()*runTime.deltaTValue();
    meanCoNum = (sum(amaxSfbyDelta)/sum(mesh.magSf())).value() * runTime.deltaTValue();
}
Foam::Info << "Mean and max Courant Numbers = "  << meanCoNum << " " << CoNum << Foam::endl;
