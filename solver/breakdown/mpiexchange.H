template<typename TScalar, typename TVector>
inline void MPIexchange
(
    const Field<TScalar> & RhoInternal,
    const Field<TVector> & UInternal,
    const Field<TScalar> & PInternal,
    const Field<TScalar> & EInternal,
    const Field<TScalar> & MolWeightInternal,
    const Field<TScalar> & InductionInternal,
    const Field<TScalar> & gammaInternal,
    const Field<TScalar> & cInternal,
    const Field<TScalar> & deltaInternal,
    Field<TScalar> & RhoNeighbour,
    Field<TVector> & UNeighbour,
    Field<TScalar> & PNeighbour,
    Field<TScalar> & ENeighbour,
    Field<TScalar> & MolWeightNeighbour,
    Field<TScalar> & InductionNeighbour,
    Field<TScalar> & gammaNeighbour,
    Field<TScalar> & cNeighbour,
    Field<TScalar> & deltaNeighbour,
    //----------------------------------------------------------------------------------------------
    const Foam::processorPolyPatch & ProcessorPolyPatch
)
{
    // делаем resize и заполняем мусором
    RhoNeighbour = RhoInternal;
    UNeighbour = UInternal;
    PNeighbour = PInternal;
    ENeighbour = EInternal;
    MolWeightNeighbour = MolWeightInternal;
    InductionNeighbour = InductionInternal;
    gammaNeighbour = gammaInternal;
    cNeighbour = cInternal;
    deltaNeighbour = deltaInternal;
    if(Pstream::master(ProcessorPolyPatch.comm()))
    {
        // Плотность
        UIPstream::read
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<char*>(RhoNeighbour.begin()),
            RhoNeighbour.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        // Скорость
        UIPstream::read
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<char*>(UNeighbour.begin()),
            UNeighbour.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        // Давление
        UIPstream::read
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<char*>(PNeighbour.begin()),
            PNeighbour.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        // Энергия
        UIPstream::read
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<char*>(ENeighbour.begin()),
            ENeighbour.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        // Молярная масса
        UIPstream::read
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<char*>(MolWeightNeighbour.begin()),
            MolWeightNeighbour.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        // Индукция
        UIPstream::read
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<char*>(InductionNeighbour.begin()),
            InductionNeighbour.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        // Показатель политропы гамма
        UIPstream::read
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<char*>(gammaNeighbour.begin()),
            gammaNeighbour.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        // Скорость звука
        UIPstream::read
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<char*>(cNeighbour.begin()),
            cNeighbour.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        // размер ячейки
        UIPstream::read
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<char*>(deltaNeighbour.begin()),
            deltaNeighbour.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        //----------------------------------------------------------------------------------------------
        // Плотность
        UOPstream::write
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<const char*>(RhoInternal.begin()),
            RhoInternal.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        // Скорость
        UOPstream::write
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<const char*>(UInternal.begin()),
            UInternal.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        // Давление
        UOPstream::write
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<const char*>(PInternal.begin()),
            PInternal.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        // Энергия
        UOPstream::write
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<const char*>(EInternal.begin()),
            EInternal.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        // Молярная масса
        UOPstream::write
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<const char*>(MolWeightInternal.begin()),
            MolWeightInternal.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        // Индукция
        UOPstream::write
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<const char*>(InductionInternal.begin()),
            InductionInternal.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        // Показатель политропы гамма
        UOPstream::write
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<const char*>(gammaInternal.begin()),
            gammaInternal.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        // Скорость звука
        UOPstream::write
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<const char*>(cInternal.begin()),
            cInternal.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        // размер ячейки
        UOPstream::write
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<const char*>(deltaInternal.begin()),
            deltaInternal.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
    }else //----------------------------------------------------------------------------------------------
    {
        // Плотность
        UOPstream::write
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<const char*>(RhoInternal.begin()),
            RhoInternal.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        // Скорость
        UOPstream::write
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<const char*>(UInternal.begin()),
            UInternal.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        // Давление
        UOPstream::write
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<const char*>(PInternal.begin()),
            PInternal.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        // Энергия
        UOPstream::write
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<const char*>(EInternal.begin()),
            EInternal.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        // Молярная масса
        UOPstream::write
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<const char*>(MolWeightInternal.begin()),
            MolWeightInternal.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        // Индукция
        UOPstream::write
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<const char*>(InductionInternal.begin()),
            InductionInternal.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        // Показатель политропы гамма
        UOPstream::write
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<const char*>(gammaInternal.begin()),
            gammaInternal.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        // Скорость звука
        UOPstream::write
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<const char*>(cInternal.begin()),
            cInternal.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        // размер ячейки
        UOPstream::write
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<const char*>(deltaInternal.begin()),
            deltaInternal.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        //----------------------------------------------------------------------------------------------
        // Плотность
        UIPstream::read
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<char*>(RhoNeighbour.begin()),
            RhoNeighbour.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        // Скорость
        UIPstream::read
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<char*>(UNeighbour.begin()),
            UNeighbour.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        // Давление
        UIPstream::read
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<char*>(PNeighbour.begin()),
            PNeighbour.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        // Энергия
        UIPstream::read
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<char*>(ENeighbour.begin()),
            ENeighbour.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        // Молярная масса
        UIPstream::read
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<char*>(MolWeightNeighbour.begin()),
            MolWeightNeighbour.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        // Индукция
        UIPstream::read
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<char*>(InductionNeighbour.begin()),
            InductionNeighbour.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        // Показатель политропы гамма
        UIPstream::read
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<char*>(gammaNeighbour.begin()),
            gammaNeighbour.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        // Скорость звука
        UIPstream::read
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<char*>(cNeighbour.begin()),
            cNeighbour.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
        // размер ячейки
        UIPstream::read
        (
            Pstream::commsTypes::blocking,
            ProcessorPolyPatch.neighbProcNo(),
            reinterpret_cast<char*>(deltaNeighbour.begin()),
            deltaNeighbour.byteSize(),
            ProcessorPolyPatch.tag(),
            ProcessorPolyPatch.comm()
        );
    }
};