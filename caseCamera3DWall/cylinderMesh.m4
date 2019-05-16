/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  4.x                                   |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

changecom(//)changequote([,])
define(calc, [esyscmd(perl -e 'printf ($1)')])
define(VCOUNT, 0)
define(vlabel, [[// ]Vertex $1 = VCOUNT define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])

meshGenApp blockMesh;
convertToMeters 0.001;


// Настройки
define(Rtr, 1.0)// Коэфициент перед периодом
define(delta, 23)// Толщина зазора
define(diameter, 306)// Внешний диаметр
//----------------------------------------------------------
define(Rlarge, calc(Rtr * diameter / 2))// Максимальный радиус
define(Rsmall, calc(Rlarge - delta))// Минимальный радиус
define(Length, 665)// Длина цилиндра
define(NPerimeter, 200)// Количество ячееу по периметру
define(NRadius, 5)// Количество ячееу вдоль радиуса
define(NLength, 120)// Количество ячееу вдоль оси
// Конец настроек

define(PI, 3.1415926535897)
define(PositionXsmall, calc(Rsmall*cos(PI/4)))
define(PositionYsmall, calc(Rsmall*sin(PI/4)))
define(PositionXlarge, calc(Rlarge*cos(PI/4)))
define(PositionYlarge, calc(Rlarge*sin(PI/4)))
define(NSection, calc(NPerimeter/4))
/*
blockMesh параметры:
Минимальный радиус - Rsmall мм
Максимальный радиус - Rlarge мм
Длина цилиндра - Length мм
Количество ячееу по периметру - NPerimeter
Количество ячееу вдоль радиуса - NRadius
Количество ячееу вдоль оси - NLength
*/
vertices
(
    ( Rsmall  0.0 0.0) // 0
    ( 0.0  Rsmall 0.0) // 1
    (-Rsmall  0.0 0.0) // 2
    ( 0.0 -Rsmall 0.0) // 3
    ( Rsmall  0.0 0.0) // 4

    ( Rlarge  0.0 0.0) // 5
    ( 0.0  Rlarge 0.0) // 6
    (-Rlarge  0.0 0.0) // 7
    ( 0.0 -Rlarge 0.0) // 8
    ( Rlarge  0.0 0.0) // 9

    ( Rsmall  0.0 Length) // 10
    ( 0.0  Rsmall Length) // 11
    (-Rsmall  0.0 Length) // 12
    ( 0.0 -Rsmall Length) // 13
    ( Rsmall  0.0 Length) // 14

    ( Rlarge  0.0 Length) // 15
    ( 0.0  Rlarge Length) // 16
    (-Rlarge  0.0 Length) // 17
    ( 0.0 -Rlarge Length) // 18
    ( Rlarge  0.0 Length) // 19
);
blocks
(
    // I четверть - верх-право
    hex (0 5 6 1 10 15 16 11) blockUPRIGHT (NRadius NSection NLength) simpleGrading (1 1 1)
    // II четверть - верх-лево
    hex (1 6 7 2 11 16 17 12) blockUPLEFT (NRadius NSection NLength) simpleGrading (1 1 1)
    // III четверть - низ-лево
    hex (2 7 8 3 12 17 18 13) blockDOWNLEFT (NRadius NSection NLength) simpleGrading (1 1 1)
    // IV четверть - низ-право
    hex (3 8 9 4 13 18 19 14) blockDOWNRIGHT (NRadius NSection NLength) simpleGrading (1 1 1)
);

//create the quarter circles
edges
(
    // Дуги для Z = 0
    // I четверть - верх-право
    arc 0 1 ( PositionXsmall  PositionYsmall 0.0)
    arc 5 6 ( PositionXlarge  PositionYlarge 0.0)
    // II четверть - верх-лево
    arc 1 2 (-PositionXsmall  PositionYsmall 0.0)
    arc 6 7 (-PositionXlarge  PositionYlarge 0.0)
    // III четверть - низ-лево
    arc 2 3 (-PositionXsmall -PositionYsmall 0.0)
    arc 7 8 (-PositionXlarge -PositionYlarge 0.0)
    // IV четверть - низ-право
    arc 3 4 ( PositionXsmall -PositionYsmall 0.0)
    arc 8 9 ( PositionXlarge -PositionYlarge 0.0)

    // Дуги для Z = Length
    // I четверть - верх-право
    arc 10 11  ( PositionXsmall  PositionYsmall Length)
    arc 15 16( PositionXlarge  PositionYlarge Length)
    // II четверть - верх-лево
    arc 11 12 (-PositionXsmall  PositionYsmall Length)
    arc 16 17(-PositionXlarge  PositionYlarge Length)
    // III четверть - низ-лево
    arc 12 13(-PositionXsmall -PositionYsmall Length)
    arc 17 18(-PositionXlarge -PositionYlarge Length)
    // IV четверть - низ-право
    arc 13 14 ( PositionXsmall -PositionYsmall Length)
    arc 18 19( PositionXlarge -PositionYlarge Length)
);

boundary
(
    inlet // Когда Z = Length
    {
        type patch;
        faces
        (
            (10 15 16 11)
            (11 16 17 12)
            (12 17 18 13)
            (13 18 19 14)
        );
    }
    outlet // Когда Z = 0
    {
        type patch;
        faces
        (
            (0 5 6 1)
            (1 6 7 2)
            (2 7 8 3)
            (3 8 9 4)
        );
    }
    walls
    {
        type wall;
        faces
        (
            // I четверть - верх-право
            (0 1 11 10)
            (5 6 16 15)
            // II четверть - верх-лево
            (1 2 12 11)
            (6 7 17 16)
            // III четверть - низ-лево
            (2 3 13 12)
            (7 8 18 17)
            // IV четверть - низ-право
            (3 4 14 13)
            (8 9 19 18)
        );
    }
    timewallUP
    {
        type patch;
        faces
        (
            (0 5 15 10)
        );
    }
    timewallDOWN
    {
        type patch;
        faces
        (
            (4 9 19 14)
        );
    }
);
mergePatchPairs
(
);