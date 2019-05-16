#define END \
{\
	if(Inversion)\
	{\
		U = -U;\
		W = -W;\
		index = -index;\
	};\
	if(ABS(U) < Epsilon*P) U = 0.0;\
	if(Ro < 0.0) Info << "Rho < 0 in face!" << endl;\
	return index;\
}

//--------------------------------------------------------------------------------------
template<class TScalar> inline void FF(TScalar P, TScalar p, TScalar ro, TScalar gamma, TScalar & F, TScalar & Fderivative)
{
	register TScalar c = SQRT(gamma*p/ro);
	if(P >= p)
	{
		F = (P-p) / (ro*c*SQRT(ABS((gamma+1.0)/gamma/2.0*P/p + (gamma-1.0)/2.0/gamma)));
		Fderivative = ((gamma + TScalar(1.0))*P/p + (TScalar(3.0)*gamma - TScalar(1.0))) / (TScalar(4.0)*gamma*ro*c*POWER((gamma+1.0)*0.5/gamma*P/p + (gamma-1.0)/gamma/TScalar(2.0), TScalar(3.0/2.0)));
	}else
	{
		F = 2.0*c/(gamma-1.0)*(POWER(ABS(P/p),(gamma-1.0)*0.5/gamma) - 1.0);
		Fderivative = c/gamma/P*POWER(ABS(P/p),(gamma-1.0)/2.0/gamma);
	}
}
//--------------------------------------------------------------------------------------
template<class TScalar> inline TScalar F(TScalar P, TScalar p, TScalar ro, TScalar gamma)
{
	TScalar c = SQRT(gamma*p/ro);
	if(P >= p) return (P-p) / (ro*c*SQRT((gamma+1.0)/gamma/2.0*P/p + (gamma-1.0)/gamma/2.0));
	else return 2.0/(gamma-1.0)*c*(POWER(ABS(P/p),(gamma-1.0)/gamma/2.0) - 1.0);
};
//--------------------------------------------------------------------------------------
template<class TScalar> inline TScalar Fderivative(TScalar P, TScalar p, TScalar ro, TScalar gamma)
{
	TScalar c = SQRT(gamma*p/ro);
	if(P >= p) return ((gamma + TScalar(1.0))*P/p + (TScalar(3.0) * gamma - TScalar(1.0))) / (TScalar(4.0)*gamma*ro*c*POWER((gamma + TScalar(1.0))/gamma/TScalar(2.0)*P/p + (gamma-1.0)/gamma/TScalar(2.0),TScalar(3.0/2.0)));
	else return c/gamma/P*POWER(ABS(P/p),(gamma-1.0)/TScalar(2.0)/gamma);
};
//--------------------------------------------------------------------------------------
template<typename TScalar>
inline int SLN(
	TScalar p1, TScalar ro1, TScalar u1, TScalar gamma1, TScalar c1,
	TScalar p2, TScalar ro2, TScalar u2, TScalar gamma2, TScalar c2,
	TScalar & P,TScalar & Ro,TScalar & U,TScalar & Ustar,TScalar W)
{
	P = 0.0;
	Ro = 0.0;
	U = 0.0;
	Ustar = 0.0;
	if(ro1 < 0.0) return 0;
	if(ro2 < 0.0) return 0;
	if(p1 < 0.0) return 0;
	if(p2 < 0.0) return 0;
	//P1 <= P2!
	bool Inversion = false;
	if(p1 > p2)
	{
		SWAP(p1,p2);
		SWAP(ro1,ro2);
		SWAP(u1,u2);
		SWAP(gamma1,gamma2);
		SWAP(c1,c2);
		u1 = -u1;
		u2 = -u2;
		W = -W;
		Inversion = true;
	};
	int index = 0;
	bool Shock1, Shock2;
	TScalar dU = u1 - u2;//Разность скоростей
	register TScalar Ushock = (p2 - p1) / SQRT(ro1*0.5*((gamma1+1.0)*p2 + (gamma1-1.0)*p1));
	TScalar Urarefaction = 2.0*c2/(gamma2-1.0)*(POWER(p1/p2,(gamma2-1.0)/gamma2/2.0) - 1.0);
	TScalar Uvacuum = -2.0*c1/(gamma1-1.0) -2.0*c2/(gamma2-1.0);

	// Определим одну из 5 возможных конфигураций и запишем в переменную index
	// 1 - Вправо и влево распространяются две ударные волны
	// 2 - Вправо - ударная волна, влево - волна разрежения
	// 3 - Влево - ударная волна, вправо - волна разрежения
	// 4 - Возникают две волны разрежения
	// 5 - Возникает область вакуума
	// P.S Случай №2 исключен, так как взяли условие, что P1 <= P2, поэтому он не рассматривается
	if(dU >= Ushock) index = 1;
	else if((Urarefaction < dU)&&(dU < Ushock)) index = 3;
	else if((Uvacuum < dU)&&(dU <= Urarefaction)) index = 4;
	else if(dU <= Uvacuum) index = 5;
	else index = 0;
	switch(index)
	{
		case 1 : Shock1 = Shock2 = true;		break;
		case 2 : Shock1 = false;Shock2 = true;	break;
		case 3 : Shock1 = true;Shock2 = false;	break;
		case 4 : Shock1 = Shock2 = false;		break;
		case 5 : index = (u1 > W) ? -index : +index; P = 1.0; Ro = 0.1; U = 0.0; END;
		default : return 0;
	}
	// Для варианта №4 давление ищется отдельно
	if(index == 4) P = p1*POWER((dU - Uvacuum)/(Urarefaction - Uvacuum),  2.0*gamma1/(gamma1-1.0));
	else
	{
		// Итерактивное вычисление давления по методу Ньютона(см. Годунова стр. 112)
		// Начальное значение заносим в P
		// Замечние: начальное приближение для исомого корня дается "снизу", т.е. не меньше
		P = (p1*ro2*c2 + p2*ro1*c1 + dU*ro1*c1*ro2*c2) / (ro1*c1 + ro2*c2);
		register TScalar Pi; //вспомогательная итерационная величина - ячейка временного значения давления
		if(P < 0.0) P = 0.0; else
		do
		{
			Pi = P;
			P = Pi - (F(Pi,p1,ro1,gamma1)+F(Pi,p2,ro2,gamma2) - dU) / (Fderivative(Pi,p1,ro1,gamma1)+Fderivative(Pi,p2,ro2,gamma2));
		}while(ABS(Pi-P) > Epsilon*Pi);
	};
	if(P < 0.0) P = 0.0;	// Давление "снизу" не может быть меньше нуля
	// После итерационного вычисления находим скорость контактного разрыва
	// Для начала вычислим массовые скорости, где a = ro(U-D), где D - это скорость ударной волны
	register TScalar a1 = (P >= p1) ? SQRT(ro1*((gamma1+1.0)/2.0*P + (gamma1-1.0)/2.0*p1)) : (ABS(p1 - P) > Epsilon*P) ? (gamma1-1.0)/2.0/gamma1*ro1*c1*(1.0-P/p1) / (1.0 - POWER(P/p1,(gamma1-1.0)/2.0/gamma1)) : ro1*c1*POWER(P/p1,(gamma1+1.0)/2.0/gamma1);
	register TScalar a2 = (P >= p2) ? SQRT(ro2*((gamma2+1.0)/2.0*P + (gamma2-1.0)/2.0*p2)) : (ABS(p2 - P) > Epsilon*P) ? (gamma2-1.0)/2.0/gamma2*ro2*c2*(1.0-P/p2) / (1.0 - POWER(P/p2,(gamma2-1.0)/2.0/gamma2)) : ro2*c2*POWER(P/p2,(gamma2+1.0)/2.0/gamma2);
	U = Ustar = (a1*u1 + a2*u2 + p1-p2) / (a1+a2);	// скорость контактного разрыва
	TScalar c1star = c1 + (gamma1 - 1.0)/2.0*(u1-U);
	TScalar c2star = c2 - (gamma2 - 1.0)/2.0*(u2-U);
	register TScalar D1 = (Shock1) ? u1 - a1/ro1 : u1 - c1;	// D1 - ударная волна влево
	register TScalar D2 = (Shock2) ? u2 + a2/ro2 : c2 + u2;	// D2 - ударная волна вправо
	register TScalar Ro1 = (Shock1) ? ro1*a1 / (a1 - ro1*(u1 - U)) : gamma1*P/(c1star*c1star);	// Плотность по адиабате Гюгонио Ro1
	register TScalar Ro2 = (Shock2) ? ro2*a2 / (a2 + ro2*(u2 - U)) : gamma2*P/(c2star*c2star);	// Плотность по адиабате Гюгонио Ro2
	TScalar D1star = U - c1star;
	TScalar D2star = U + c2star;
	//------------------------------------------------------
	// Надо найти, в какую из зон попадает величина W-движение грани
	//_______________________________________________________________________________________________________________> W
	//                   |                    |                     |                    |                    |
	//                   D1                  D1star                 U                   D2star                D2  
	
	if(W <= D1)								// 1 случай по Годунову на стр. 120 - берем всю левую область(u1,p1,ro1)
	{
		Ro = ro1;
		P = p1;
		U = u1;
		index = -index;
		END;
	}else if ((D1 < W) && (W < D1star))		// 4(1) случай - продолжение 1 случая по Годунову на стр. 120
	{
		register TScalar cstar = (gamma1 - 1.0) / (gamma1 + 1.0)*(u1 - W) + 2.0/(gamma1 + 1.0)*c1;
		U = W + cstar;
		P = p1*POWER(ABS(cstar)/c1, 2.0*gamma1 / (gamma1 - 1.0));
		Ro = gamma1*P / (cstar*cstar);
		index = -index;
		END;
	}
	else if(W >= D2)						// 2 случай по Годунову на стр. 120 - берем всю правую область(u2,p2,ro2)
	{
		Ro = ro2;
		P = p2;
		U = u2;
		index = +index;
		END;
	}else if ((D2star < W) && (W < D2))		// 4(2) случай - продолжение 2 случая по Годунову на стр. 120
	{
		register TScalar cstar = (gamma2 - 1.0) / (gamma2 + 1.0)*(u2 - W) - 2.0/(gamma2 + 1.0)*c2;
		U = W + cstar;
		P = p2*POWER(ABS(cstar)/c2, 2.0*gamma2 / (gamma2 - 1.0));
		Ro = gamma2*P / (cstar*cstar);
		index = +index;
		END;
	}
	else if((D1star <= W)&&(W <= U))		// 3 случай по Годунову на стр. 120 - берем область(U,P,Ro1)
	{
		Ro = Ro1;
		// U,P уже заданы
		index = -index;
		END;
	}else if((U <= W) && (W <= D2star))		// 3 случай по Годунову на стр. 120 - берем область(U,P,Ro2)
	{
		Ro = Ro2;
		// U,P уже заданы
		index = +index;
		END;
	}else
	{
		return 0;	// Ошибка
	}
    END;
};
//=================================================
template<typename TScalar>
inline int SLN_accustic(
	TScalar p1, TScalar ro1, TScalar u1, TScalar gamma1, TScalar c1,
	TScalar p2, TScalar ro2, TScalar u2, TScalar gamma2, TScalar c2,
	TScalar & P,TScalar & Ro,TScalar & U,TScalar & Ustar,TScalar W)
{
	int index = 0;
	P = (p1*ro2*c2 + p2*ro1*c1 + (u1-u2)*ro1*c1*ro2*c2) / (ro1*c1 + ro2*c2);
	U = u1 - (P-p1)/(ro1*c1);
	Ustar = U;
	if(Ustar > W)
	{
		Ro = ro1 + (P-p1)/(c1*c1);
		index = -1;
	}else
	{
		Ro = ro2 + (P-p2)/(c2*c2);
		index = 1;
	}
    return index;
};
