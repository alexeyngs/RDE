// Список скалярных переменных
// 0 - время сущствования жесткой стенки
// 1 - SstarS
// 2 - Pstar
// 3 - Tstar
// 4 - Backbressure
Foam::scalarIOList scalarParameters
(
    IOobject
    (
        "scalarParameters",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    5
);

Foam::dictionary dictionaryBoundary("./constant/Properties");
std::filebuf fb;
bool b = fb.open(dictionaryBoundary.name(), std::ios::in);
std::istream stream(&fb);
Foam::ISstream Sstream(stream, dictionaryBoundary.name());
dictionaryBoundary.read(Sstream);
// Чтение
scalarParameters[0] = 0.0;
scalarParameters[1] = readScalar(dictionaryBoundary.subDict("geometry").lookup("SstarS"));
scalarParameters[2] = readScalar(dictionaryBoundary.subDict("geometry").lookup("Pstar"));
scalarParameters[3] = readScalar(dictionaryBoundary.subDict("geometry").lookup("Tstar"));
scalarParameters[4] = readScalar(dictionaryBoundary.subDict("geometry").lookup("Backbressure"));
// Закрытие
fb.close();
