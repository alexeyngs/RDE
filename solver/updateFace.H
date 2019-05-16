gamma = thermo.Cp()/thermo.Cv();
// Скорость звука
c = sqrt(gamma * rPsi);

// Потоки уравнений неразрывности
//---------Start---------
RhoUf = (Rhof*Uf) & mesh.Sf();
//---------End---------

// Потоки уравнений импульса
//---------Start---------
RhoUUf = (Rhof*Uf*Uf) & mesh.Sf();
PPf = Pf * mesh.Sf();
//---------End---------

// Потоки уравнений энергии
//---------Start---------
Hf = Ef + Pf/Rhof + 0.5*magSqr(Uf);
RhoUHf = (Rhof*Uf*Hf) & mesh.Sf();
//---------End---------

// Потоки уравнений молярной массы
//---------Start---------
RhoUMolWeightf = (Rhof*Uf*MolWeightf) & mesh.Sf();
//---------End---------

// Потоки уравнений доли индукции
//---------Start---------
RhoUInductionf = (Rhof*Uf*Inductionf) & mesh.Sf();
//---------End---------


// ребра при линейной интерполяции
Rholf = linearInterpolate(rho);
Ulf = linearInterpolate(U);
Plf = linearInterpolate(p);
Tlf = linearInterpolate(T);
Elf = linearInterpolate(e);
MolWeightlf = linearInterpolate(MolWeight);
Inductionlf = linearInterpolate(Induction);
clf = linearInterpolate(c);
gammalf = linearInterpolate(gamma);
