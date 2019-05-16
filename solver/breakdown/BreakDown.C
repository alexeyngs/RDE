#include "fvCFD.H"
#include "tmp.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "surfaceInterpolationScheme.H"
#include "basicThermo.H"
#include "pointPatch.H"
#include "processorFvPatch.H"
#include "../thermo/rhoRDEThermo.H"
template<class TScalar> inline TScalar ABS(TScalar X)
{
	return X > 0 ? X : -X;
};
template<class TScalar> TScalar SQRT(TScalar X);
template<class TScalar> TScalar POWER(TScalar X, TScalar Y); // X^Y
template<class TScalar> inline void SWAP(TScalar & X, TScalar & Y)
{
	TScalar buffer = X;
	X = Y;
	Y = buffer;
};
extern bool TimeWall;
const Foam::scalar Epsilon = Foam::SMALL*10.0;
#include "sln.H"
#include "mpiexchange.H"
#include "Godunov.H"
#include "Kolgan.H"
#include "Accustic.H"
enum class EBreakDown
{
    None,
    Accustic,
    Godunov,
    Kolgan,
};
EBreakDown TypeBreakDown = EBreakDown::None;
void ReadBreakDown(const Foam::dictionary & dict)
{
    Foam::word BreakDownString = dict.lookup("BreakDown");
    if(BreakDownString == "Accustic") TypeBreakDown = EBreakDown::Accustic;
    else if(BreakDownString == "Godunov") TypeBreakDown = EBreakDown::Godunov;
    else if(BreakDownString == "Kolgan") TypeBreakDown = EBreakDown::Kolgan;
}
void BreakDown(
    Foam::rhoRDEThermo & Thermo, const Foam::volScalarField & Psi,
    const Foam::volScalarField & Rho, const Foam::volVectorField & U, const Foam::volScalarField & P, const Foam::volScalarField & E, const Foam::volScalarField & MolWeight, const Foam::volScalarField & Induction,
    const Foam::surfaceScalarField & Rholf, const Foam::surfaceVectorField & Ulf, const Foam::surfaceScalarField & Plf, const Foam::surfaceScalarField & Elf, const Foam::surfaceScalarField & MolWeightlf, const Foam::surfaceScalarField & Inductionlf, const Foam::surfaceScalarField & clf, const Foam::surfaceScalarField & gammalf,
    Foam::surfaceScalarField & Rhof, Foam::surfaceVectorField & Uf, Foam::surfaceScalarField & Pf, Foam::surfaceScalarField & Ef, Foam::surfaceScalarField & Ethermodynamicalf, Foam::surfaceScalarField & MolWeightf, Foam::surfaceScalarField & Inductionf)
{
    Foam::volScalarField gamma = Thermo.Cp() / Thermo.Cv();
//    volScalarField c = sqrt(gamma / Psi);
    Foam::volScalarField c = sqrt(gamma * P / Rho);
    switch(TypeBreakDown)
    {
        case EBreakDown::None : Foam::Perr << endl << "NO BREAKDOWN!"<< Foam::endl; abort(); return;
        case EBreakDown::Accustic : Accustic(Thermo, Psi, geometricOneField(), Rho, U, P, E, MolWeight, Induction, c, gamma,  Rhof, Uf, Pf, Ef, Ethermodynamicalf, MolWeightf, Inductionf);break;
        case EBreakDown::Godunov :   Godunov(Thermo, Psi, geometricOneField(), Rho, U, P, E, MolWeight, Induction, c, gamma, Rhof, Uf, Pf, Ef, Ethermodynamicalf, MolWeightf, Inductionf);break;
        case EBreakDown::Kolgan :     Kolgan(Thermo, Psi, geometricOneField(), Rho, U, P, E, MolWeight, Induction, c, gamma, Rholf, Ulf, Plf, Elf, MolWeightlf, Inductionlf, clf, gammalf, Rhof, Uf, Pf, Ef, Ethermodynamicalf, MolWeightf, Inductionf);break;
        default : return;
    };
}










// Переделать...
//         ----------------------------------------------------------
//         |             |             ||             |             |
//         |     L1x     |     L1      ||      L2     |     L2x     |
//         |             |             ||             |             |
//         ----------------------------------------------------------
void BreakDown(
    const Foam::vector Normal, const Foam::scalar L1, const Foam::scalar L1x, const Foam::scalar L2, const Foam::scalar L2x,
    const Foam::scalar Rho1, const Foam::vector U1, const Foam::scalar P1, const Foam::scalar E1, const Foam::scalar MolWeight1, const Foam::scalar Induction1, const Foam::scalar gamma1,
    const Foam::scalar Rho2, const Foam::vector U2, const Foam::scalar P2, const Foam::scalar E2, const Foam::scalar MolWeight2, const Foam::scalar Induction2, const Foam::scalar gamma2,
    const Foam::scalar Rho1x, const Foam::vector U1x, const Foam::scalar P1x, const Foam::scalar E1x, const Foam::scalar MolWeight1x, const Foam::scalar Induction1x, const Foam::scalar gamma1x,
    const Foam::scalar Rho2x, const Foam::vector U2x, const Foam::scalar P2x, const Foam::scalar E2x, const Foam::scalar MolWeight2x, const Foam::scalar Induction2x, const Foam::scalar gamma2x,
    Foam::scalar & Rhof, Foam::vector & Uf, Foam::scalar & Pf, Foam::scalar & Ef, Foam::scalar & MolWeightf, Foam::scalar & Inductionf)
{
    Foam::scalar c1 = Foam::sqrt(gamma1 * P1 / Rho1);
    Foam::scalar c2 = Foam::sqrt(gamma2 * P2 / Rho2);
    switch(TypeBreakDown)
    {
        case EBreakDown::None : Foam::Perr << endl << "NO BREAKDOWN!"<< Foam::endl; abort(); return;
        case EBreakDown::Accustic : 
        {
            // Разложение скорости
            Foam::scalar UU1 = U1 & Normal;         // Проекция скорости на нормаль - нормальная компонента скорости в ячейке 1
		    Foam::vector T1 = U1 - (Normal * UU1);	// Тангенсальная компонента скорости в ячейке 1        
            Foam::scalar UU2 = U2 & Normal;		    // Проекция скорости на нормаль - нормальная компонента скорости в ячейке 2
            Foam::vector T2 = U2 - (Normal * UU2);	// Тангенсальная компонента скорости в ячейке 2
            //---------------------------------------------
            Foam::scalar Ustar;                      // скорость контактного разрыва
            Foam::scalar W = 0.0;
            Foam::scalar UUf;
            int sln = SLN_accustic<Foam::scalar>(P1, Rho1, UU1, gamma1, c1, P2, Rho2, UU2, gamma2, c2, Pf, Rhof, UUf, Ustar, W);
            if(sln == -5 || sln == +5) Info << "Vacuum" << endl;
            if(sln < 0)
            {
                Uf = T1 + (Normal * UUf);
                MolWeightf = MolWeight1;
                Inductionf = Induction1;
                Ef = E1;
            }else if(sln > 0)
            {
                Uf = T2 + (Normal * UUf);
                MolWeightf = MolWeight2;
                Inductionf = Induction2;
                Ef = E2;
            }else
            {
                Info << "ERORR!" << endl;
                break;	//Error
            }
        }break;

        case EBreakDown::Godunov :
        {
            // Разложение скорости
            Foam::scalar UU1 = U1 & Normal;         // Проекция скорости на нормаль - нормальная компонента скорости в ячейке 1
		    Foam::vector T1 = U1 - (Normal * UU1);	// Тангенсальная компонента скорости в ячейке 1        
            Foam::scalar UU2 = U2 & Normal;		    // Проекция скорости на нормаль - нормальная компонента скорости в ячейке 2
            Foam::vector T2 = U2 - (Normal * UU2);	// Тангенсальная компонента скорости в ячейке 2
            //---------------------------------------------
            Foam::scalar Ustar;                      // скорость контактного разрыва
            Foam::scalar W = 0.0;
            Foam::scalar UUf;
            int sln = SLN<Foam::scalar>(P1, Rho1, UU1, gamma1, c1, P2, Rho2, UU2, gamma2, c2, Pf, Rhof, UUf, Ustar, W);
            if(sln == -5 || sln == +5) Info << "Vacuum" << endl;
            if(sln < 0)
            {
                Uf = T1 + (Normal * UUf);
                MolWeightf = MolWeight1;
                Inductionf = Induction1;
                Ef = E1;
            }else if(sln > 0)
            {
                Uf = T2 + (Normal * UUf);
                MolWeightf = MolWeight2;
                Inductionf = Induction2;
                Ef = E2;
            }else
            {
                Info << "ERORR!" << endl;
                break;	//Error
            }
        }break;

        case EBreakDown::Kolgan :
        {
            Foam::scalar P1_ = P1 + gradient(L1x, L1, P1x, (P1x*L1 + P1*L1x) / (L1x+L1), P1);
            Foam::scalar Rho1_ = Rho1 + gradient(L1x, L1, Rho1x, (Rho1x*L1 + Rho1*L1x) / (L1x+L1), Rho1);
            Foam::vector U1_ = U1 + gradient(L1x, L1, U1x, (U1x*L1 + U1*L1x) / (L1x+L1), U1);
            Foam::scalar gamma1_ = gamma1 + gradient(L1x, L1, gamma1x, (gamma1x*L1 + gamma1*L1x) / (L1x+L1), gamma1);
            Foam::scalar E1_ = E1 + gradient(L1x, L1, E1x, (E1x*L1 + E1*L1x) / (L1x+L1), E1);
            Foam::scalar MolWeight1_ = MolWeight1 + gradient(L1x, L1, MolWeight1x, (MolWeight1x*L1 + MolWeight1*L1x) / (L1x+L1), MolWeight1);
            Foam::scalar Induction1_ = Induction1 + gradient(L1x, L1, Induction1x, (Induction1x*L1 + Induction1*L1x) / (L1x+L1), Induction1);

            Foam::scalar P2_ = P2 - gradient(L2, L2x, P2, (P2*L2x + P2x*L2) / (L2 + L2x), P2x);
            Foam::scalar Rho2_ = Rho2 - gradient(L2, L2x, Rho2, (Rho2*L2x + Rho2x*L2) / (L2 + L2x), Rho2x);
            Foam::vector U2_ = U2 - gradient(L2, L2x, U2, (U2*L2x + U2x*L2) / (L2 + L2x), U2x);
            Foam::scalar gamma2_ = gamma2 - gradient(L2, L2x, gamma2, (gamma2*L2x + gamma2x*L2) / (L2 + L2x), gamma2x);
            Foam::scalar E2_ = E2 - gradient(L2, L2x, E2, (E2*L2x + E2x*L2) / (L2 + L2x), E2x);
            Foam::scalar MolWeight2_ = MolWeight2 - gradient(L2, L2x, MolWeight2, (MolWeight2*L2x + MolWeight2x*L2) / (L2 + L2x), MolWeight2x);
            Foam::scalar Induction2_ = P2 - gradient(L2, L2x, Induction2, (Induction2*L2x + Induction2x*L2) / (L2 + L2x), Induction2x);

            // Разложение скорости
            Foam::scalar UU1 = U1_ & Normal;        // Проекция скорости на нормаль - нормальная компонента скорости в ячейке 1
		    Foam::vector T1 = U1_ - (Normal * UU1);	// Тангенсальная компонента скорости в ячейке 1        
            Foam::scalar UU2 = U2_ & Normal;        // Проекция скорости на нормаль - нормальная компонента скорости в ячейке 2
            Foam::vector T2 = U2_ - (Normal * UU2);	// Тангенсальная компонента скорости в ячейке 2
            //---------------------------------------------
            Foam::scalar Ustar;                     // скорость контактного разрыва
            Foam::scalar W = 0.0;
            Foam::scalar UUf;
            int sln = SLN<Foam::scalar>(P1_, Rho1_, UU1, gamma1_, c1, P2_, Rho2_, UU2, gamma2_, c2, Pf, Rhof, UUf, Ustar, W);
            if(sln == -5 || sln == +5) Info << "Vacuum" << endl;
            if(sln < 0)
            {
                Uf = T1 + (Normal * UUf);
                MolWeightf = MolWeight1_;
                Inductionf = Induction1_;
                Ef = E1_;
            }else if(sln > 0)
            {
                Uf = T2 + (Normal * UUf);
                MolWeightf = MolWeight2_;
                Inductionf = Induction2_;
                Ef = E2_;
            }else
            {
                Info << "ERORR!" << endl;
                break;	//Error
            }
        }break;
        default : return;
    };
}