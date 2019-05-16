
#ifndef InletDictionary_H
#define InletDictionary_H

#include "scalarIOList.H"

Foam::scalar Pin(const Foam::scalarIOList & scalarParameters, const Foam::scalar time, const Foam::scalar molWeight, const Foam::scalar Pcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar gamma);
Foam::scalar Tin(const Foam::scalarIOList & scalarParameters, const Foam::scalar time, const Foam::scalar molWeight, const Foam::scalar Pcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar gamma);
Foam::vector Uin(const Foam::scalarIOList & scalarParameters, const Foam::scalar time, const Foam::scalar molWeight, const Foam::scalar Pcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar gamma);
Foam::scalar Rhoin(const Foam::scalarIOList & scalarParameters, const Foam::scalar time, const Foam::scalar molWeight, const Foam::scalar Pcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar gamma);

Foam::scalar Pout(const Foam::scalarIOList & scalarParameters, const Foam::scalar time, const Foam::scalar molWeight, const Foam::scalar Pcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar gamma);
Foam::scalar Tout(const Foam::scalarIOList & scalarParameters, const Foam::scalar time, const Foam::scalar molWeight, const Foam::scalar Pcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar gamma);
Foam::vector Uout(const Foam::scalarIOList & scalarParameters, const Foam::scalar time, const Foam::scalar molWeight, const Foam::scalar Pcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar gamma);
Foam::scalar Rhoout(const Foam::scalarIOList & scalarParameters, const Foam::scalar time, const Foam::scalar molWeight, const Foam::scalar Pcenter, const Foam::scalar Rhocenter, const Foam::vector Ucenter, const Foam::scalar gamma);

#endif