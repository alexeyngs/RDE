#include "makeReactionThermo.H"
#include "rhoReactionThermo.H"
#include "heRhoThermo.H"
#include "specie.H"
#include "perfectGas.H"
#include "incompressiblePerfectGas.H"
#include "thermo.H"

#include "constTransport.H"
#include "sutherlandTransport.H"

#include "absoluteEnthalpy.H"
#include "sensibleEnthalpy.H"


#include "RDEspecie.H"
#include "rhoRDEThermo.H"
#include "heRDEThermo.H"
#include "RDEThermo.H"


#include "RDEMixtureBasic.H"
#include "RDEMixtureNikolaev.H"
#include "RDEMixtureVasilev.H"
#include "RDEMixtureEquilibrium.H"
namespace Foam
{
    // без химии
    makeReactionThermo
    (
        rhoRDEThermo,
        heRDEThermo,
        RDEMixtureNikolaev,
        constTransport,
        sensibleEnthalpy,
        RDEThermo,
        perfectGas,
        RDEspecie
    );
    // Николаев
    makeReactionThermo
    (
        rhoRDEThermo,
        heRDEThermo,
        RDEMixtureNikolaev,
        constTransport,
        absoluteEnthalpy,
        RDEThermo,
        perfectGas,
        RDEspecie
    );
    // Васильев
    makeReactionThermo
    (
        rhoRDEThermo,
        heRDEThermo,
        RDEMixtureVasilev,
        constTransport,
        absoluteEnthalpy,
        RDEThermo,
        perfectGas,
        RDEspecie
    );
    // Равновесная темодинамика
    makeReactionThermo
    (
        rhoRDEThermo,
        heRDEThermo,
        RDEMixtureEquilibrium,
        constTransport,
        absoluteEnthalpy,
        RDEThermo,
        perfectGas,
        RDEspecie
    );
}
// Foam::species::thermo<Thermo, Type> = Foam::species::thermo<Foam::RDEThermo<Foam::perfectGas<Foam::RDEspecie> >, Foam::absoluteEnthalpy>
// В mixture: ThermoType = Foam::constTransport<Foam::species::thermo<Foam::RDEThermo<Foam::perfectGas<Foam::RDEspecie> >, Foam::sensibleEnthalpy> >
// Thermo = Foam::RDEThermo<Foam::perfectGas<Foam::RDEspecie> >;
// Type = Foam::absoluteEnthalpy;
// Foam::scalar = double;
// Foam::label = int


/*
Код OpenFOAM очень широко использует шаблоны. При этом очень активно используется следующая интересная структура:
Есть, допустим, некоторый базовый класс. И есть целое множество шаблонных классов, наследующих свой параметр:
template <class T> class Foo: public T

при этом каждый такой класс добавляет к своему предку какие-либо новые свойства/особенности,
а в итоге используется целая цепочка таких шаблонов типа Foo<Bar<Baz > >.

Например, простейший класс, хранящий термодинамические свойства вещества, имеет вид constTransport<specieThermo<hConstThermo > >.
Здесь perfectGas — класс, отвечающий (грубо) за вычисление плотности газа по давлению и температуре и т.п.;
hConstThermo добавляет к классу хранение (постоянной) теплоемкости при фиксированном давлении и вычисление энтальпии по температуре и наоборот;
specieThermo позволяет пересчитывать другие функции типа теплоемкости при фиксированном объеме или энтропии;
constTransport добавляет хранение (постоянных) вязкости и теплопроводности и предоставляет к ним доступ.

Для каждого такого «уровня» есть, как правило, несколько реализаций. Например, вместо constTransport можно использовать sutherlandTransport,
вычисляющий вязкость и теплопроводность по температуре в соответствии с формулой Сазерленда.
В результате можно сочетать разные модели для разных термодинамических свойств, тем самым порождая те комбинации, которые требуются.
Аналогичные конструкции встречаются во всем коде OpenFOAM.
*/


/*
const Foam::dimensionSet Foam::dimless(0, 0, 0, 0, 0, 0, 0);

const Foam::dimensionSet Foam::dimMass(1, 0, 0, 0, 0, 0, 0);
const Foam::dimensionSet Foam::dimLength(0, 1, 0, 0, 0, 0, 0);
const Foam::dimensionSet Foam::dimTime(0, 0, 1, 0, 0, 0, 0);
const Foam::dimensionSet Foam::dimTemperature(0, 0, 0, 1, 0, 0, 0);
const Foam::dimensionSet Foam::dimMoles(0, 0, 0, 0, 1, 0, 0);
const Foam::dimensionSet Foam::dimCurrent(0, 0, 0, 0, 0, 1, 0);
const Foam::dimensionSet Foam::dimLuminousIntensity(0, 0, 0, 0, 0, 0, 1);

const Foam::dimensionSet Foam::dimArea(sqr(dimLength));
const Foam::dimensionSet Foam::dimVolume(pow3(dimLength));
const Foam::dimensionSet Foam::dimVol(dimVolume);

const Foam::dimensionSet Foam::dimVelocity(dimLength/dimTime);
const Foam::dimensionSet Foam::dimAcceleration(dimVelocity/dimTime);

const Foam::dimensionSet Foam::dimDensity(dimMass/dimVolume);
const Foam::dimensionSet Foam::dimForce(dimMass*dimAcceleration);
const Foam::dimensionSet Foam::dimEnergy(dimForce*dimLength);
const Foam::dimensionSet Foam::dimPower(dimEnergy/dimTime);

const Foam::dimensionSet Foam::dimPressure(dimForce/dimArea);
const Foam::dimensionSet Foam::dimGasConstant(dimEnergy/dimMass/dimTemperature);
const Foam::dimensionSet Foam::dimSpecificHeatCapacity(dimGasConstant);
*/
