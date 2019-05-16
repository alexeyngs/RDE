template<class ThermoType>
inline Foam::scalar GetMolWeightInduction(const speciesTable & species, const PtrList<ThermoType> & speciesData, const Foam::scalar Induction, const Foam::scalar powerN, const Foam::scalar Sigma)
{
    const Foam::scalar BETA = Foam::pow(Induction, powerN);
    const scalar molesCH4 = speciesData[species["CH4"]].moles();
    const scalar molesH2 = speciesData[species["H2"]].moles();
    const scalar molesO2 = speciesData[species["O2"]].moles();
    const scalar molesN2 = speciesData[species["N2"]].moles();
    const scalar molesAr = speciesData[species["Ar"]].moles();
    const scalar molesH2O = speciesData[species["H2O"]].moles();
    //--------------------------------------------------------------
    Foam::scalar MolWeight = Sigma / (2.5*molesCH4 - 1.5*BETA*molesCH4 + molesH2 + molesO2 + molesN2 + molesAr + molesH2O);
    return MolWeight;
}

template<class ThermoType>
inline Foam::scalar GetHchemical(const speciesTable & species, const PtrList<ThermoType> & speciesData, const Foam::scalar Induction, const Foam::scalar powerN, const Foam::scalar Sigma)
{
    const Foam::scalar BETA = Foam::pow(Induction, powerN);
    const scalar molesCH4 = speciesData[species["CH4"]].moles();
    const scalar molesH2 = speciesData[species["H2"]].moles();
    const scalar molesO2 = speciesData[species["O2"]].moles();
    const scalar molesN2 = speciesData[species["N2"]].moles();
    const scalar molesH2O = speciesData[species["H2O"]].moles();
    //--------------------------------------------------------------
    const Foam::scalar Ed_CH4 = speciesData[species["CH4"]].Edissociation_;
    const Foam::scalar Ed_H2 = speciesData[species["H2"]].Edissociation_;
    const Foam::scalar Ed_O2 = speciesData[species["O2"]].Edissociation_;
    const Foam::scalar Ed_H2O = speciesData[species["H2O"]].Edissociation_;
    const Foam::scalar Ed_OH = speciesData[species["OH"]].Edissociation_;
    Foam::scalar Hchemical = (BETA*molesCH4*Ed_CH4 + (1.0 - BETA)*molesCH4*Ed_OH + (1.5*molesCH4 - 1.5*BETA*molesCH4 + molesH2)*Ed_H2 + (molesO2 - molesCH4 + BETA*molesCH4)*Ed_O2 + molesH2O*Ed_H2O) / Sigma;
    return Hchemical*1000.0;
}

template<class ThermoType>
inline Foam::scalar GetMolWeightMax(const speciesTable & species, const PtrList<ThermoType> & speciesData, const Foam::scalar fuel, const Foam::scalar Sigma)
{
    // Количество молей в смеси
    const Foam::scalar molesCH4 = speciesData[species["CH4"]].moles();
    const Foam::scalar molesH2 = speciesData[species["H2"]].moles();
    const Foam::scalar molesO2 = speciesData[species["O2"]].moles();
    const Foam::scalar molesN2 = speciesData[species["N2"]].moles();
    const Foam::scalar molesAr = speciesData[species["Ar"]].moles();
    const Foam::scalar molesH2O = speciesData[species["H2O"]].moles();
    // Количество атомов в смеси
    const Foam::scalar nC = molesCH4;
    const Foam::scalar nH = 4.0*molesCH4 + 2.0*molesH2 + 2.0*molesH2O;
    const Foam::scalar nO = 2.0*molesO2 + molesH2O;
    const Foam::scalar nN = 2.0*molesN2;
    const Foam::scalar nAr = molesAr;
    //--------------------------------------------------------------
    Foam::scalar molWeightMax = 0.0;
    if(nC <= nO < 2.0*nC) molWeightMax = Sigma / (nC + nH/2.0 + nN/2.0 + nAr);
    if(2.0*nC <= nO < 2.0*nH) molWeightMax = Sigma / (nC + nH/2.0 + nN/2.0 + nAr);
    if(nO >= 2.0*nC + nH/2.0) molWeightMax = Sigma / (nH/4.0 + nO/2.0 + nN/2.0 + nAr);
    return molWeightMax;
}

template<class ThermoType>
inline Foam::scalar GetEdissociation(const speciesTable & species, const PtrList<ThermoType> & speciesData, const Foam::scalar fuel, const Foam::scalar Sigma, const Foam::scalar molWeightMin)
{
    const Foam::scalar molesCH4 = speciesData[species["CH4"]].moles();
    const Foam::scalar molesH2 = speciesData[species["H2"]].moles();
    const Foam::scalar molesO2 = speciesData[species["O2"]].moles();
    const Foam::scalar molesN2 = speciesData[species["N2"]].moles();
    const Foam::scalar molesAr = speciesData[species["Ar"]].moles();
    const Foam::scalar molesH2O = speciesData[species["H2O"]].moles();
    //--------------------------------------------------------------
    const Foam::scalar Ed_CH4 = speciesData[species["CH4"]].Edissociation_;
    const Foam::scalar Ed_H2 = speciesData[species["H2"]].Edissociation_;
    const Foam::scalar Ed_O2 = speciesData[species["O2"]].Edissociation_;
    const Foam::scalar Ed_H2O = speciesData[species["H2O"]].Edissociation_;
    const Foam::scalar Ed_CO2 = speciesData[species["CO2"]].Edissociation_;
    const Foam::scalar Ed_OH = speciesData[species["OH"]].Edissociation_;
    const Foam::scalar Ed_CO = 0.0; // Точка отсчета от CO, поэтому EdCO равен 0
    Foam::scalar Qmax = 0.0;
    Foam::scalar Qstar = 0.0;
    //--------------------------------------------------------------
    if(molesCH4/2.0 <= molesO2 && molesO2 < molesCH4)
    {
        // Кислорода хватило только на окисление углерода и немного водорода
        // Разложение топлива CH4 + O2 -> CO + H2
        Qstar = (2.0*molesO2 - molesCH4)*Ed_CO + (2.0*molesO2 - molesCH4)*Ed_OH + (2.5*molesCH4 + molesH2 - molesO2)*Ed_H2;
    }else
    {
        // Кислорода хватит полностью на CO и ОН
        // Разложение топлива CH4 + O2 -> CO + OH + 1.5H2
        Qstar = molesCH4*Ed_CO + molesCH4*Ed_OH + 1.5*molesCH4*Ed_H2 + (molesO2 - molesCH4)*Ed_O2 + molesH2*Ed_H2;
    }
    //--------------------------------------------------------------
    if(molesO2 >= 2.0*molesCH4 + molesH2/2.0)   // Явно избыток кислорода
    {
        Qmax = molesCH4*Ed_CO2 + (2.0*molesCH4)*Ed_H2O + molesH2*Ed_H2O + (molesO2 - 2.0*molesCH4 - molesH2/2.0)*Ed_O2;
    }else if(molesCH4 <= molesO2 && molesO2 < 2.0*molesCH4 + molesH2/2.0)   // Кислорода хватит только на окисления углерода CO2
    {
        Qmax = molesCH4*Ed_CO2 + (2.0*molesO2 - 2.0*molesCH4)*Ed_H2O + (4.0*molesCH4 + molesH2 - 2.0*molesO2)*Ed_H2;
    }else if(molesO2 < molesCH4)    // Кислорода не хватило даже на полное окисление углерода, но на образование CO кислорода есть всегда, т.к. должны получать газ
    {
        Qmax = (2.0*molesO2 - molesCH4)*Ed_CO2 + 2.0*molesCH4*Ed_H2;
    };
    register Foam::scalar Q0 = molesCH4*Ed_CH4 + molesO2*Ed_O2 + molesH2*Ed_H2; // Отсчет от mu0
    Qmax = -Qmax;
    Q0 = -Q0;
    Qstar = -Qstar;
    Info << "molesO2 = " << molesO2 << ", Qmax = " << Qmax/4.1868 << "(" << (Qmax-Q0)/4.1868<< ")" << " cal, Q* = " << Qstar/4.1868 << "(" << (Qstar-Q0)/4.1868<< ")" << "cal , Q0 = " << Q0/4.1868 << " cal" << endl;
    Foam::scalar molWeightStar = GetMolWeightInduction(species, speciesData, 0.0, 1.0, Sigma);
    Foam::scalar molWeightMax = GetMolWeightMax(species, speciesData, fuel, Sigma);
    // Энергия диссоциации от * до max:
//    Foam::scalar Edissociation = (Qmax - Qstar) / Sigma / (1.0/molWeightStar - 1.0/molWeightMax);
    // Энергия диссоциации от min до max:
    Foam::scalar Edissociation = Qmax / Sigma / (1.0/molWeightMin - 1.0/molWeightMax);
    return Edissociation;
}


/*
    if(fuel > 2.0*molesO2)  // Если обогащенная смесь
    {
    }else   // Если обедненная смесь
    {
    }
    // Текущая энергия диссоциации для конкретной ячейки
    scalar Edi = 0.0;
    forAll(speciesData, i)
    {
        Edi += speciesData[i].Edissociation_ * speciesData[i].moles();
    }
*/