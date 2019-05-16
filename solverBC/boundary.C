#include <iostream>     // std::ios, std::istream, std::cout
#include <fstream>
#include "fixedValueFvPatchFields.H"
#include "thermodynamicConstants.H"
#include "fvCFD.H"
#include "scalarIOList.H"

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
    word BreakDownString = dict.lookup("BreakDown");
    if(BreakDownString == "Accustic") TypeBreakDown = EBreakDown::Accustic;
    else if(BreakDownString == "Godunov") TypeBreakDown = EBreakDown::Godunov;
    else if(BreakDownString == "Kolgan") TypeBreakDown = EBreakDown::Kolgan;
}
//====================================================================================================================================

// Входящие параметры
Foam::scalar Pin(const Foam::scalarIOList & scalarParameters, const scalar time, const Foam::scalar molWeight, const Foam::scalar Pcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar gamma)
{
    const Foam::scalar & timeWall = scalarParameters[0];
    const Foam::scalar & SstarS = scalarParameters[1];
    const Foam::scalar & Pstar = scalarParameters[2];
    const Foam::scalar & Tstar = scalarParameters[3];
    const Foam::scalar & Backbressure = scalarParameters[4];
    //------------------------------------------------------------------------
    return Pcenter;
}

Foam::vector Uin(const Foam::scalarIOList & scalarParameters, const scalar time, const Foam::scalar molWeight, const Foam::scalar Pcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar gamma)
{
    const Foam::scalar & timeWall = scalarParameters[0];
    const Foam::scalar & SstarS = scalarParameters[1];
    const Foam::scalar & Pstar = scalarParameters[2];
    const Foam::scalar & Tstar = scalarParameters[3];
    const Foam::scalar & Backbressure = scalarParameters[4];
    //------------------------------------------------------------------------
    Foam::scalar F, X, U;
    if(Pcenter >= Pstar)
    {
        U = 0.0;
        return Foam::vector(Ucenter.x(), Ucenter.y(), 0.0);
    }else
    {
        F = SstarS * Foam::sqrt(gamma * Foam::pow(2.0/(gamma + 1.0), (gamma + 1.0)/(gamma - 1.0)));
        X = gamma/(gamma - 1.0) * Pcenter/Pstar/F;
        U = X - Foam::sqrt(X*X + 2.0*gamma/(gamma - 1.0));
        //----------------------
        Foam::scalar Rhostar = Pstar * molWeight / Foam::constant::thermodynamic::RR / Tstar;
        U *= Foam::sqrt(Pstar/Rhostar); // Выход из нормировки на *
	}
    return Foam::vector(0.0, 0.0, U);
}

Foam::scalar Rhoin(const Foam::scalarIOList & scalarParameters, const scalar time, const Foam::scalar molWeight, const Foam::scalar Pcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar gamma)
{
    const Foam::scalar & timeWall = scalarParameters[0];
    const Foam::scalar & SstarS = scalarParameters[1];
    const Foam::scalar & Pstar = scalarParameters[2];
    const Foam::scalar & Tstar = scalarParameters[3];
    const Foam::scalar & Backbressure = scalarParameters[4];
    //------------------------------------------------------------------------
    Foam::scalar rho;
    Foam::scalar F, X, U;
    if(Pcenter >= Pstar)
    {
        rho = Rhocenter;
    }else
    {
        F = SstarS * Foam::sqrt(gamma * Foam::pow(2.0/(gamma + 1.0), (gamma + 1.0)/(gamma - 1.0)));
        X = gamma/(gamma - 1.0) * Pcenter/Pstar/F;
        U = X - Foam::sqrt(X*X + 2.0*gamma/(gamma - 1.0));
        //----------------------
        Foam::scalar Rhostar = Pstar * molWeight / Foam::constant::thermodynamic::RR / Tstar;
        rho = -F/U;
        rho *= Rhostar; // Выход из нормировки на *
	}
    return rho;
}

Foam::scalar Tin(const Foam::scalarIOList & scalarParameters, const scalar time, const Foam::scalar molWeight, const Foam::scalar Pcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar gamma)
{
    const Foam::scalar & timeWall = scalarParameters[0];
    const Foam::scalar & SstarS = scalarParameters[1];
    const Foam::scalar & Pstar = scalarParameters[2];
    const Foam::scalar & Tstar = scalarParameters[3];
    const Foam::scalar & Backbressure = scalarParameters[4];
    //------------------------------------------------------------------------
    Foam::scalar Rho = Rhoin(scalarParameters, time, molWeight, Pcenter, Rhocenter, Ucenter, gamma);
    Foam::scalar P = Pin(scalarParameters, time, molWeight, Pcenter, Rhocenter, Ucenter, gamma);
    return P * molWeight / Rho / Foam::constant::thermodynamic::RR;
}

// Выходящие параметры
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
const Foam::scalar Epsilon = Foam::SMALL*10.0;
#include "../solver/breakdown/sln.H"
Foam::scalar Pout(const Foam::scalarIOList & scalarParameters, const scalar time, const Foam::scalar molWeight, const Foam::scalar Pcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar gamma)
{
    const Foam::scalar & timeWall = scalarParameters[0];
    const Foam::scalar & SstarS = scalarParameters[1];
    const Foam::scalar & Pstar = scalarParameters[2];
    const Foam::scalar & Tstar = scalarParameters[3];
    const Foam::scalar & Backbressure = scalarParameters[4];
    //------------------------------------------------------------------------
    if(time < timeWall) return Pcenter;
    Foam::scalar W = 0.0;
    Foam::scalar P, Rho, U, Ustar;
    int result = SLN<Foam::scalar>(Backbressure, Rhocenter, Ucenter.z(), gamma, Foam::sqrt(gamma*Backbressure/Rhocenter), Pcenter, Rhocenter, Ucenter.z(), gamma, Foam::sqrt(gamma*Pcenter/Rhocenter), P, Rho, U, Ustar, W);
    return P;
}

Foam::scalar Tout(const Foam::scalarIOList & scalarParameters, const scalar time, const Foam::scalar molWeight, const Foam::scalar Pcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar gamma)
{
    const Foam::scalar & timeWall = scalarParameters[0];
    const Foam::scalar & SstarS = scalarParameters[1];
    const Foam::scalar & Pstar = scalarParameters[2];
    const Foam::scalar & Tstar = scalarParameters[3];
    const Foam::scalar & Backbressure = scalarParameters[4];
    //------------------------------------------------------------------------
    if(time < timeWall) return Pcenter * molWeight / Rhocenter / Foam::constant::thermodynamic::RR;
    Foam::scalar W = 0.0;
    Foam::scalar P, Rho, U, Ustar;
    int result = SLN<Foam::scalar>(Backbressure, Rhocenter, Ucenter.z(), gamma, Foam::sqrt(gamma*Backbressure/Rhocenter), Pcenter, Rhocenter, Ucenter.z(), gamma, Foam::sqrt(gamma*Pcenter/Rhocenter), P, Rho, U, Ustar, W);
    return P * molWeight / Rho / Foam::constant::thermodynamic::RR;
}

Foam::vector Uout(const Foam::scalarIOList & scalarParameters, const scalar time, const Foam::scalar molWeight, const Foam::scalar Pcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar gamma)
{
    const Foam::scalar & timeWall = scalarParameters[0];
    const Foam::scalar & SstarS = scalarParameters[1];
    const Foam::scalar & Pstar = scalarParameters[2];
    const Foam::scalar & Tstar = scalarParameters[3];
    const Foam::scalar & Backbressure = scalarParameters[4];
    //------------------------------------------------------------------------
    if(time < timeWall) return Ucenter;
    Foam::scalar W = 0.0;
    Foam::scalar P, Rho, U, Ustar;
    int result = SLN<Foam::scalar>(Backbressure, Rhocenter, Ucenter.z(), gamma, Foam::sqrt(gamma*Backbressure/Rhocenter), Pcenter, Rhocenter, Ucenter.z(), gamma, Foam::sqrt(gamma*Pcenter/Rhocenter), P, Rho, U, Ustar, W);
    return Foam::vector(Ucenter.x(), Ucenter.y(), U);
}

Foam::scalar Rhoout(const Foam::scalarIOList & scalarParameters, const scalar time, const Foam::scalar molWeight, const Foam::scalar Pcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar gamma)
{
    const Foam::scalar & timeWall = scalarParameters[0];
    const Foam::scalar & SstarS = scalarParameters[1];
    const Foam::scalar & Pstar = scalarParameters[2];
    const Foam::scalar & Tstar = scalarParameters[3];
    const Foam::scalar & Backbressure = scalarParameters[4];
    //------------------------------------------------------------------------
    if(time < timeWall) return Rhocenter;
    Foam::scalar W = 0.0;
    Foam::scalar P, Rho, U, Ustar;
    int result = SLN<Foam::scalar>(Backbressure, Rhocenter, Ucenter.z(), gamma, Foam::sqrt(gamma*Backbressure/Rhocenter), Pcenter, Rhocenter, Ucenter.z(), gamma, Foam::sqrt(gamma*Pcenter/Rhocenter), P, Rho, U, Ustar, W);
    return Rho;
}
