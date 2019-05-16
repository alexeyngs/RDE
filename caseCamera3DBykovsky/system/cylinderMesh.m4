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

//----------------------------------------------------------
// Настройки
define(Rtr, 1.0)// Коэфициент перед периодом
define(delta, 23)// Толщина зазора камеры сгорания
define(sdelta, 2)// Добавочная толщина зазора коллектора
define(diameter, 306)// Внешний диаметр
define(CameraHeight, 665)// Длина камеры сгорания
define(CollectorHeight, 42)// Длина коллектора
//Ячейки:
define(NPerimeter, 300)// Количество ячеек по периметру
define(NRadiusCollector, 2)// Количество ячеек вдоль радиуса коллектора
define(NRadiusCamera, 23)// Количество ячеек вдоль камеры сгорания
define(NCollector, 10)// Количество ячеек вдоль оси коллектора
define(NCamera, 100)// Количество ячеек вдоль оси камеры сгорания
// Конец настроек
//----------------------------------------------------------

define(Length, calc(CameraHeight + CollectorHeight)) // Общая длина
define(Rlarge, calc(Rtr * diameter / 2)) // Максимальный радиус
define(Rmiddle, calc(Rlarge - adddelta)) // Средний радиус
define(Rsmall, calc(Rlarge - adddelta - delta)) // Минимальный радиус
define(PI, 3.1415926535897)
define(PositionXsmall, calc(Rsmall*cos(PI/4)))
define(PositionYsmall, calc(Rsmall*sin(PI/4)))
define(PositionXmiddle, calc(Rmiddle*cos(PI/4)))
define(PositionYmiddle, calc(Rmiddle*sin(PI/4)))
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
    // Вход в коллектор
    ( Rsmall  0.0 Length) // 0
    ( 0.0  Rsmall Length) // 1
    (-Rsmall  0.0 Length) // 2
    ( 0.0 -Rsmall Length) // 3

    ( Rmiddle  0.0 Length) // 4
    ( 0.0  Rmiddle Length) // 5
    (-Rmiddle  0.0 Length) // 6
    ( 0.0 -Rmiddle Length) // 7

    // Выход из коллектора - вход в камеру сгорания
    ( Rsmall  0.0 CameraHeight) // 8
    ( 0.0  Rsmall CameraHeight) // 9
    (-Rsmall  0.0 CameraHeight) // 10
    ( 0.0 -Rsmall CameraHeight) // 11

    ( Rmiddle  0.0 CameraHeight) // 12
    ( 0.0  Rmiddle CameraHeight) // 13
    (-Rmiddle  0.0 CameraHeight) // 14
    ( 0.0 -Rmiddle CameraHeight) // 15

    ( Rlarge  0.0 CameraHeight) // 16
    ( 0.0  Rlarge CameraHeight) // 17
    (-Rlarge  0.0 CameraHeight) // 18
    ( 0.0 -Rlarge CameraHeight) // 19

    // Выход из камеры сгорания
    ( Rsmall  0.0 0.0) // 20
    ( 0.0  Rsmall 0.0) // 21
    (-Rsmall  0.0 0.0) // 22
    ( 0.0 -Rsmall 0.0) // 23

    ( Rmiddle  0.0 0.0) // 24
    ( 0.0  Rmiddle 0.0) // 25
    (-Rmiddle  0.0 0.0) // 26
    ( 0.0 -Rmiddle 0.0) // 27

    ( Rlarge  0.0 0.0) // 28
    ( 0.0  Rlarge 0.0) // 29
    (-Rlarge  0.0 0.0) // 30
    ( 0.0 -Rlarge 0.0) // 31
);

blocks
(
    // Коллектор
    // I четверть - верх-право
    hex (8 12 13 9 0 4 5 1) blockUPRIGHT (NCollector NSection NRadiusCollector) simpleGrading (1 1 1)
    // II четверть - верх-лево
    hex (9 13 14 10 1 5 6 2) blockUPLEFT (NCollector NSection NRadiusCollector) simpleGrading (1 1 1)
    // III четверть - низ-лево
    hex (10 14 15 11 2 6 7 3) blockDOWNLEFT (NCollector NSection NRadiusCollector) simpleGrading (1 1 1)
    // IV четверть - низ-право
    hex (11 15 12 8 3 7 4 0) blockDOWNRIGHT (NCollector NSection NRadiusCollector) simpleGrading (1 1 1)

    // Камера сгорания
    // I четверть - верх-право
//    hex (8 16 17 9 20 28 29 21) blockUPRIGHT (NCamera NSection NRadiusCamera) simpleGrading (1 1 1)
    // II четверть - верх-лево
//    hex (9 17 18 10 21 29 30 22) blockUPLEFT (NCamera NSection NRadiusCamera) simpleGrading (1 1 1)
    // III четверть - низ-лево
//    hex (10 18 19 11 22 30 31 23) blockDOWNLEFT (NCamera NSection NRadiusCamera) simpleGrading (1 1 1)
    // IV четверть - низ-право
//    hex (11 19 16 8 23 31 28 20) blockDOWNRIGHT (NCamera NSection NRadiusCamera) simpleGrading (1 1 1)
);

//create the quarter circles
edges
(
    // Дуги для Z = вход в коллектор
    // I четверть - верх-право
    arc 0 1 ( PositionXsmall  PositionYsmall 0.0)
    arc 4 5 ( PositionXmiddle  PositionYmiddle 0.0)
    // II четверть - верх-лево
    arc 1 2 (-PositionXsmall  PositionYsmall 0.0)
    arc 5 6 (-PositionXmiddle  PositionYmiddle 0.0)
    // III четверть - низ-лево
    arc 2 3 (-PositionXsmall -PositionYsmall 0.0)
    arc 6 7 (-PositionXmiddle -PositionYmiddle 0.0)
    // IV четверть - низ-право
    arc 3 0 ( PositionXsmall -PositionYsmall 0.0)
    arc 7 4 ( PositionXmiddle -PositionYmiddle 0.0)

    // Дуги для Z = вход в камеру сгорания
    // I четверть - верх-право
    arc 8 9  ( PositionXsmall  PositionYsmall CameraHeight)
    arc 12 13( PositionXmiddle  PositionYmiddle CameraHeight)
//    arc 16 17( PositionXlarge  PositionYlarge CameraHeight)
    // II четверть - верх-лево
    arc 9 10 (-PositionXsmall  PositionYsmall CameraHeight)
    arc 13 14(-PositionXmiddle  PositionYmiddle CameraHeight)
//    arc 17 18(-PositionXlarge  PositionYlarge CameraHeight)
    // III четверть - низ-лево
    arc 10 11(-PositionXsmall -PositionYsmall CameraHeight)
    arc 14 15(-PositionXmiddle -PositionYmiddle CameraHeight)
//    arc 18 19(-PositionXlarge -PositionYlarge CameraHeight)
    // IV четверть - низ-право
    arc 11 8 ( PositionXsmall -PositionYsmall CameraHeight)
    arc 15 12( PositionXmiddle -PositionYmiddle CameraHeight)
//    arc 19 16( PositionXlarge -PositionYlarge CameraHeight)
);

boundary
(
    inlet // Когда Z = Length
    {
        type patch;
        faces
        (
            (0 4 5 1)
            (1 5 6 2)
            (2 6 7 3)
            (3 7 4 0)
        );
    }
    outlet // Когда Z = 0
    {
        type patch;
        faces
        (
            (8 12 13 9)
            (9 13 14 10)
            (10 14 15 11)
            (11 15 12 8)
        );
    }
    walls
    {
        type wall;
        faces
        (
            // I четверть - верх-право
            (0 1 9 8)
            (4 5 13 12)
            // II четверть - верх-лево
            (1 2 10 9)
            (5 6 14 13)
            // III четверть - низ-лево
            (2 3 11 10)
            (6 7 15 14)
            // IV четверть - низ-право
            (3 0 8 11)
            (7 4 12 15)
        );
    }
/*
    walls
    {
        type wall;
        faces
        (
            // I четверть - верх-право
            (12 16 17 13)
            // II четверть - верх-лево
            (13 17 18 14)
            // III четверть - низ-лево
            (14 18 19 15)
            // IV четверть - низ-право
            (15 19 16 12)
        );
    }
*/
);
mergePatchPairs
(
);