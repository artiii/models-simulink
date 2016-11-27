model TransferFunction "Linear transfer function" 
  parameter Real b[:]={1} 
    "Numerator coefficients of transfer function.";
  parameter Real a[:]={1,1} 
    "Denominator coefficients of transfer function.";
  output Real x[size(a, 1) - 1] 
    "State of transfer function from controller canonical form";
protected 
  parameter Integer na=size(a, 1) 
    "Size of Denominator of transfer function.";
  parameter Integer nb(max=na) = size(b, 1) 
    "Size of Numerator of transfer function.";
  parameter Integer nx=size(a, 1) - 1;
  Real x1dot "Derivative of first state of TransferFcn";
  Real xn "Highest order state of TransferFcn";
  Real u; 
  Real y; 
equation 
  [der(x); xn] = [x1dot; x];
  [u] = transpose([a])*[x1dot; x];
  [y] = transpose([zeros(na - nb, 1); b])*[x1dot; x];
end TransferFunction;

model ExcitationRegulator "Simple model of ARV-SDP1"
  parameter Real K0u=-10;
  parameter Real K1u=0;
  parameter Real K0w=0;
  parameter Real K1w=0;
  parameter Real K1if=0;
  parameter Real DL0=0;
  constant Real Tokp=0.05;
  constant Real T0u=0.02;
  constant Real T1u=0.039;
  constant Real Tfb=0.07;
  constant Real T0w=1.0;
  constant Real T1w=0.026;
  constant Real T1if=0.03;
  TransferFunction VoltageDeviation(b={K0u}, a={T0u, 1});
  TransferFunction VoltageDerivative(b={K1u, 0}, a={T1u, 1});
  TransferFunction FrequencyBlockD(b={1},   a={Tfb, 1});
  TransferFunction FrequencyBlockU(b={1,0}, a={Tfb, 1});
  TransferFunction FrequencyDeviation(b={K0w, 0}, a={T0w, 1});
  TransferFunction FrequencyDerivative(b={K1w, 0}, a={T1w, 1});
  TransferFunction ExcitationCurrent(b={K1if, 0}, a={T1if, 1});
  TransferFunction SimpleExciter(b={1}, a={Tokp, 1});
  Real u_u, u_pu, u_fsys, u_fu, u_pIf, y_er, f_sum;
equation
  when initial() then
     reinit(FrequencyBlockU.x[1],DL0);
//     reinit(FrequencyDeviation.x[1],DL0);
//     reinit(FrequencyDerivative.x[1],DL0);
  end when;
  u_u = VoltageDeviation.u;
  u_pu = VoltageDerivative.u;
  u_fsys = FrequencyBlockD.u;
  u_fu   = FrequencyBlockU.u;
  u_pIf = ExcitationCurrent.u;
  f_sum = FrequencyBlockD.y + FrequencyBlockU.y;
//  connect(FrequencyBlock.y, FrequencyDeviation.u);
//  connect(FrequencyBlock.y, FrequencyDerivative.u);
  FrequencyDeviation.u  = f_sum;
  FrequencyDerivative.u = f_sum;
  VoltageDeviation.y + VoltageDerivative.y + ExcitationCurrent.y +
      FrequencyDeviation.y + FrequencyDerivative.y = SimpleExciter.u;
  y_er = SimpleExciter.y;
end ExcitationRegulator;

model Constant_Conductivity_Load "The representation of load as a constant conductivity"
//-------------------------------------------------------------
  parameter Real Gn=0.7 "Активная составляющая нагрузки";
  parameter Real Bn=-0.35 "Реактивная составляющая нагрузки";
//-------------------------------------------------------------
  parameter Real TloadOff=1000;
  parameter Real dTloadOff=10;
//-------------------------------------------------------------
  Real Udn, Uqn, Iqn, Idn;
  Real Pn, Qn;
  Real Ssys;
  Real R0, X0, Gnn, Bnn;
//  constant Real pi=4*arctan(1);
//  constant Real wc=100*pi;
equation
  when initial() then
     R0 =  Gn/(Gn^2+Bn^2);
     X0 = -Bn/(Gn^2+Bn^2);
  end when;
  Gnn = Gn; //            R0/(R0^2 + (X0*(1+Ssys))^2);
  Bnn = Bn; //-(X0*(1+Ssys))/(R0^2 + (X0*(1+Ssys))^2);
  Pn =  Uqn*Iqn + Udn*Idn;
  Qn = -Uqn*Idn + Udn*Iqn;
  Idn = if time>=TloadOff and time<TloadOff+dTloadOff then 0
                                                      else Udn*Gnn + Uqn*Bnn;
  Iqn = if time>=TloadOff and time<TloadOff+dTloadOff then 0
                                                      else Uqn*Gnn - Udn*Bnn;
end Constant_Conductivity_Load;

model ShortCircuitShunt "Three phase short circuit analogue"
//-----------------------------------
  parameter Real Bn=-100000;
//-----------------------------------
  parameter Real TkzOn=1000;
  parameter Real dTkzOn=0.12;
//-----------------------------------
  Real id, iq, Ud, Uq;
equation
  id = if time>=TkzOn and time<TkzOn+dTkzOn then  Bn*Uq
                                            else 0;
  iq = if time>=TkzOn and time<TkzOn+dTkzOn then -Bn*Ud
                                            else 0;
end ShortCircuitShunt;

model Transformer "The model of HV transformer"
//---------------------------------------------
  parameter Real Rline=0.001;
  parameter Real Xline=0.01;
//---------------------------------------------
  parameter Real TtOff=1000;
  parameter Real dTtOff=10;
//---------------------------------------------
  Real U1d, U1q, U2d, U2q, I12d, I12q, I12m, U1m, U2m;
  Real P1, Q1, P2, Q2, RL, XL;
equation
// the algorithm of triping transmission line
  RL = if time>=TtOff and time<TtOff+dTtOff then 100000*Rline
                                            else Rline;
  XL = if time>=TtOff and time<TtOff+dTtOff then 100000*Xline
                                            else Xline;
// series impedance between nodes
  U1q = U2q + RL*I12q - XL*I12d;
  U1d = U2d + RL*I12d + XL*I12q;
// measurements
  I12m = sqrt(I12d^2+I12q^2);
  U1m = sqrt(U1d^2+U1q^2);
  U2m = sqrt(U2d^2+U2q^2);
// powers
  P1 =  U1q*I12q + U1d*I12d;
  Q1 = -U1q*I12d + U1d*I12q;
  P2 =  U2q*I12q + U2d*I12d;
  Q2 = -U2q*I12d + U2d*I12q;
end Transformer;

model HVLine "The model of extra HV transmission line"
//---------------------------------------------
  parameter Real Rline=0.01;
  parameter Real Xline=0.1;
  parameter Real Bline=0.1;
//---------------------------------------------
  parameter Real TLineOff=1000;
  parameter Real dTLineOff=10;
//---------------------------------------------
  Real U1d, U1q, U2d, U2q, I12d, I12q, I12m, U1m, U2m;
  Real I1d, I1q, I2d, I2q, Idc1, Iqc1, Idc2, Iqc2, DU1, DU2, Ssys, Bnn;
  Real P1, Q1, P2, Q2, RL, XL;
equation
// the algorithm of triping transmission line
  RL = if time>=TLineOff and time<TLineOff+dTLineOff then 100000*Rline
                                                     else Rline;
  XL = if time>=TLineOff and time<TLineOff+dTLineOff then 100000*Xline
                                                     else Xline;
// series impedance between nodes
  U1q = U2q + RL*I12q - XL*I12d;
  U1d = U2d + RL*I12d + XL*I12q;
// line capacitances
  Bnn = Bline; // /(1+Ssys);
  Idc1 =  U1q*(Bnn/2);
  Iqc1 = -U1d*(Bnn/2);
  Idc2 =  U2q*(Bnn/2);
  Iqc2 = -U2d*(Bnn/2);
// balance of current
  I1d = Idc1 + I12d;
  I1q = Iqc1 + I12q;
  I12d = Idc2 + I2d;
  I12q = Iqc2 + I2q;
// measurements
  I12m = sqrt(I12d^2+I12q^2);
  U1m = sqrt(U1d^2+U1q^2);
  U2m = sqrt(U2d^2+U2q^2);
  DU1 = arctan2(U1d,U1q);
  DU2 = arctan2(U2d,U2q);
// powers
  P1 =  U1q*I1q + U1d*I1d;
  Q1 = -U1q*I1d + U1d*I1q;
  P2 =  U2q*I2q + U2d*I2d;
  Q2 = -U2q*I2d + U2d*I2q;
end HVLine;

model SyncronousMachine "Full model of syncronous machine" 
//------------------------------------------------------------------------------------
  parameter Real TgenOff=1000;
  parameter Real dTgenOff=10;
//------------------------------------------------------------------------------------
  parameter Real Pg=0.85;
  parameter Real Qg=0.527;
  parameter Real Ut=1;
//------------------------------------------------------------------------------------
  parameter Real Xd=1.869 "Продольное индуктивное сопротивление";
  parameter Real Xq=1.869 "Поперечное индуктивное сопротивление";
  parameter Real Xs=0.194 "Индуктивное сопротивление рассеяния";
  parameter Real X1d=0.3016 "Переходное продольное индуктивное сопротивление";
  parameter Real X2d=0.2337 "Сверхпереходное продольное индуктивное сопротивление";
  parameter Real X2q=0.2337 "Сверхпереходное поперечное индуктивное сопротивление";
  parameter Real Rf=904e-6 "Активное сопротивление обмотки возбуждения";
  parameter Real R1d=3.688e-3 "Активное сопротивление ДК в продольной оси";
  parameter Real R1q=2.77e-3 "Активное сопротивление ДК в поперечной оси";
  parameter Real Tj=7 "Механическая инерционная постоянная";
  constant  Real Xt=0 "Индуктивное сопротивление трансформатора";
  parameter Real Sigma=0.0475 "Коэффициент статизма регулятора турбины";
  parameter Real TauC=0.8 "Постоянная времени сервомотора";
//------------------------------------------------------------------------------------
  Real dUtr, dWu, dIf, dEr, dWf;
  constant Real PI=4*arctan(1);
  constant Real Wc=100*PI;
// Parameters of power system elements
  constant Real Ra=0;
  constant Real Rt=0;
  constant Real Mtmax=1.1;
  //------------------------------------------------------------------
  Real EQ, Dg, Uq, Ud, Id, Iq, Mt, Eiq, Ir, Uf;
  Real Eq, Er0, PsiR, PsiRD, PsiRQ, Er_max, Er_min, Uf_full;
// Uc, Dl, Ug, Dt, Pt, Qt, , X2ds, X2qs, Udgen, Uqgen, Ugen0,
  Real Xad, Xaq, X2dt, X2qt, Rs, Xsf, Xs1d, Xs1q, X1, X2, X3;
  Real Ugen;
// Ut, Pg, Qg, 
  Real Mu0, Ro, Mu, Mt_pp;
  Real Pgen, Qgen;
  Real UdG, UqG, IdG, IqG, DeltaIJ, Ssys, Tregmax;
  //------------------------------------------------------------------
// Integrated Variables
  Real Yr, Yrd, Yrq, s, DGi;
  Real Me, Yad, Yaq, iq, id, ir, ird, irq;
  Real ud, uq, E11d, E11q;
equation 
  when initial() then
    Xad=Xd - Xs;
    Xaq=Xq - Xs;
    X2dt=X2d + Xt;
    X2qt=X2q + Xt;
    Rs=Ra + Rt;
    Xsf=1/(1/(X1d - Xs) - 1/Xad);
    Xs1d=1/(1/(X2d - Xs) - 1/Xad - 1/Xsf);
    Xs1q=1/(1/(X2q - Xs) - 1/Xaq);
    X1=(X2d-Xs)/Xsf;
    X2=(X2d-Xs)/Xs1d;
    X3=(X2q-Xs)/Xs1q;
    EQ=sqrt((Ut + ((Ra+Rt)*Pg + (Xq+Xt)*Qg)/Ut)^2 +
           (((Xq+Xt)*Pg - (Ra+Rt)*Qg)/Ut)^2);
    Dg=arctan((((Xq+Xt)*Pg-(Ra+Rt)*Qg)/Ut)/(Ut+((Ra+Rt)*Pg+(Xq+Xt)*Qg)/Ut));
//----------------------
    Uq=Ut*cos(Dg);
    Ud=-Ut*sin(Dg);
//----------------------
    Id=((Uq-EQ)*(Xq+Xt)-(Ra+Rt)*Ud)/((Ra+Rt)^2 + (Xq+Xt)^2);
    Iq=(-Ud-(Ra+Rt)*Id)/(Xq+Xt);
    Mt=EQ*Iq;
    Mu0=Mt;
    Eiq=Uq-Id*(Xs+Xt)+(Ra+Rt)*Iq;
    Ir=-Id+Eiq/Xad;
    Uf=Ir*Rf;
    Eq=Ir*Xad;
    Er0=Eq;
    PsiR=Eiq + Xsf*Ir;
    PsiRQ=-Ud-Iq*(Xs+Xt)-Id*(Ra+Rt);
    PsiRD=Eiq;
    Er_max=2*Er0;
    Er_min=-0.6*Er0;
//-------Initialization State Variables------------------
    reinit(DGi, Dg);
    reinit(s, 0);
    reinit(Yr, PsiR);
    reinit(Yrd, PsiRD);
    reinit(Yrq, PsiRQ);
    reinit(Mu, Mu0);
  end when;
  der(DGi) = Wc*s;

  Uf_full = if      (Uf*Xad/Rf+dEr) > Er_max then Er_max*Rf/Xad
            else if (Uf*Xad/Rf+dEr) < Er_min then Er_min*Rf/Xad
                                             else Uf + Rf*dEr/Xad;

  der(Yr) = Wc*(Uf_full - Rf*ir);
  der(Yrd) = -Wc*R1d*ird;
  der(Yrq) = -Wc*R1q*irq;
  der(s) =(Mt_pp - Me)/Tj;
  Ro = -Mu + Mu0 - s/Sigma;
  der(Mu) = if time<Tregmax then Ro/TauC
                            else Ro/10e6;

  Mt_pp = if Mu>=0 and Mu<Mtmax*Pg then Mu
              else if Mu>=Mtmax*Pg then Mtmax*Pg
                                   else 0;

  Me = Yad*iq - Yaq*id;
  ir = (Yr - Yad)/Xsf;
  ird = (Yrd - Yad)/Xs1d;
  irq = (Yrq - Yaq)/Xs1q;
  id =  if time>=TgenOff and time<TgenOff+dTgenOff then 0
                         else ((uq-E11q)*X2qt-Rs*(ud+E11d))/(X2dt*X2qt+Rs^2);
  iq =  if time>=TgenOff and time<TgenOff+dTgenOff then 0
                         else (-ud-E11d-Rs*id)/X2qt;
  E11d = X3*Yrq;
  E11q = X1*Yr + X2*Yrd;
  Yad = E11q + id*(X2d - Xs);
  Yaq = E11d + iq*(X2q - Xs);
  Ugen = sqrt(ud^2+uq^2);
  dUtr = sqrt(ud^2+uq^2)-Ut;
  dWu = Wc*Ssys;
  dWf = if UqG <> 0 then arctan2(UdG,UqG)
                    else arctan2(UdG,0.001);
  dIf = ir - Ir;
//----------------------------------------
  Pgen =  uq*iq + ud*id;
  Qgen = -uq*id + ud*iq;
//-----------------------------------------------------------
  uq = UqG*cos(DeltaIJ) - UdG*sin(DeltaIJ);
  ud = UqG*sin(DeltaIJ) + UdG*cos(DeltaIJ);
//-----------------------------------------------------------
  IqG = iq*cos(DeltaIJ) + id*sin(DeltaIJ);
  IdG = id*cos(DeltaIJ) - iq*sin(DeltaIJ);
//-----------------------------------------------------------
end SyncronousMachine;

model Generator_with_ARV "The syncronous machine eqiupped with strong excitation"
  SyncronousMachine G;
  ExcitationRegulator AVR;
equation
  connect(G.dUtr,AVR.u_u);
  connect(G.dUtr,AVR.u_pu);
  connect(G.dWu, AVR.u_fsys);
  connect(G.dWf, AVR.u_fu);
  connect(G.dIf, AVR.u_pIf);
  connect(AVR.y_er,G.dEr);
//--------------proper frequency measurement------------------
//-------AVR.u_f = Delta_G1_bus + arctan2(G.ud,G.uq);---------
//------------------------------------------------------------
end Generator_with_ARV;

model Kursovik "The simplified model for term paper calculation"
  parameter Real TregMax=50;
//--------------------------------------------------------------
  Generator_with_ARV G1;
  Generator_with_ARV G2;
  Generator_with_ARV G3;

  Transformer T_G1;
  Transformer T_G2;
  Transformer T1;
  Transformer T2_HV;
  Transformer T2_MV;
  Transformer T2_LV;
  Transformer T3;
  Transformer AT_HV;
  Transformer AT_MV;
  Transformer AT_LV;

  HVLine L1;
  HVLine L2;
  HVLine L3;
  HVLine L4;
  HVLine L5;
  HVLine L6;

  Constant_Conductivity_Load N1;
  Constant_Conductivity_Load N2;
  Constant_Conductivity_Load N3;
  Constant_Conductivity_Load N4;

  ShortCircuitShunt KZ01;
  ShortCircuitShunt KZ02;
  ShortCircuitShunt KZ03;
  ShortCircuitShunt KZ04;
  ShortCircuitShunt KZ05;
  ShortCircuitShunt KZ06;
  ShortCircuitShunt KZ07;
  ShortCircuitShunt KZ08;
  ShortCircuitShunt KZ09;
  ShortCircuitShunt KZ10;
  ShortCircuitShunt KZ11;
  ShortCircuitShunt KZ12;
  ShortCircuitShunt KZ13;
  ShortCircuitShunt KZ14;
  ShortCircuitShunt KZ15;
  ShortCircuitShunt KZ16;
//--------------------------------------------
  Real delta_G1_G2, delta_G1_G3, delta_G2_G3;
equation
//-------------------------The equations of current balance----------------------------------
//--------------------------------------------(Node 01)---------------------------------------
  G2.G.IdG - T_G2.I12d - KZ01.id = 0;
  G2.G.IqG - T_G2.I12q - KZ01.iq = 0;
//--------------------------------------------(Node 02)---------------------------------------
  T_G2.I12d - L2.I1d - KZ02.id = 0;
  T_G2.I12q - L2.I1q - KZ02.iq = 0;
//--------------------------------------------(Node 03)---------------------------------------
  L2.I2d - AT_HV.I12d - KZ03.id = 0;
  L2.I2q - AT_HV.I12q - KZ03.iq = 0;
//--------------------------------------------(Node 04)---------------------------------------
  AT_HV.I12d - AT_MV.I12d - AT_LV.I12d - KZ04.id = 0;
  AT_HV.I12q - AT_MV.I12q - AT_LV.I12q - KZ04.iq = 0;
//--------------------------------------------(Node 05)---------------------------------------
  AT_MV.I12d - L6.I1d - KZ05.id = 0;
  AT_MV.I12q - L6.I1q - KZ05.iq = 0;
//--------------------------------------------(Node 06)---------------------------------------
  AT_LV.I12d - N1.Idn - KZ06.id = 0;
  AT_LV.I12q - N1.Iqn - KZ06.iq = 0;
//--------------------------------------------(Node 07)---------------------------------------
  L1.I2d + L6.I2d - L3.I1d - L4.I1d - T3.I12d - KZ07.id = 0;
  L1.I2q + L6.I2q - L3.I1q - L4.I1q - T3.I12q - KZ07.iq = 0;
//--------------------------------------------(Node 08)---------------------------------------
  T_G1.I12d - L1.I1d - KZ08.id = 0;
  T_G1.I12q - L1.I1q - KZ08.iq = 0;
//--------------------------------------------(Node 09)---------------------------------------
  G1.G.IdG - T_G1.I12d - KZ09.id = 0;
  G1.G.IqG - T_G1.I12q - KZ09.iq = 0;
//--------------------------------------------(Node 10)---------------------------------------
  L4.I2d + L5.I2d - T1.I12d - KZ10.id = 0;
  L4.I2q + L5.I2q - T1.I12q - KZ10.iq = 0;
//--------------------------------------------(Node 11)---------------------------------------
  T1.I12d - N2.Idn - KZ11.id = 0;
  T1.I12q - N2.Iqn - KZ11.iq = 0;
//--------------------------------------------(Node 12)---------------------------------------
  L3.I2d - L5.I1d - T2_MV.I12d - KZ12.id = 0;
  L3.I2q - L5.I1q - T2_MV.I12q - KZ12.iq = 0;
//--------------------------------------------(Node 13)---------------------------------------
  T2_HV.I12d + T2_MV.I12d - T2_LV.I12d - KZ13.id = 0;
  T2_HV.I12q + T2_MV.I12q - T2_LV.I12q - KZ13.iq = 0;
//--------------------------------------------(Node 14)---------------------------------------
  G3.G.IdG - T2_HV.I12d - KZ14.id = 0;
  G3.G.IqG - T2_HV.I12q - KZ14.iq = 0;
//--------------------------------------------(Node 15)---------------------------------------
  T2_LV.I12d - N3.Idn - KZ15.id = 0;
  T2_LV.I12q - N3.Iqn - KZ15.iq = 0;
//--------------------------------------------(Node 16)---------------------------------------
  T3.I12d - N4.Idn - KZ16.id = 0;
  T3.I12q - N4.Iqn - KZ16.iq = 0;
//--------------------------------------------------------------------------------------------

//-----------Angles of Machines---------------------------
  der(delta_G1_G2) = G1.G.Wc*(G1.G.s - G2.G.s);
  der(delta_G1_G3) = G1.G.Wc*(G1.G.s - G3.G.s);
  der(delta_G2_G3) = G1.G.Wc*(G2.G.s - G3.G.s);
  
  G1.G.DeltaIJ = 0;
  G2.G.DeltaIJ = delta_G1_G2;
  G3.G.DeltaIJ = delta_G1_G3;

  G1.G.Ssys=G1.G.s;
  G2.G.Ssys=G1.G.s;
  G3.G.Ssys=G1.G.s;

//---------------------------Line connections-----------------------
//--------------------------------------(Node 01)-------------------
  G2.G.UqG = KZ01.Uq;
  G2.G.UdG = KZ01.Ud;

  T_G2.U1q   = KZ01.Uq;
  T_G2.U1d   = KZ01.Ud;
//--------------------------------------(Node 02)-------------------
  T_G2.U2q = KZ02.Uq;
  T_G2.U2d = KZ02.Ud;

  L2.U1q = KZ02.Uq;
  L2.U1d = KZ02.Ud;
//--------------------------------------(Node 03)-------------------
  L2.U2q = KZ03.Uq;
  L2.U2d = KZ03.Ud;

  AT_HV.U1q = KZ03.Uq;
  AT_HV.U1d = KZ03.Ud;
//--------------------------------------(Node 04)-------------------
  AT_HV.U2q   = KZ04.Uq;
  AT_HV.U2d   = KZ04.Ud;

  AT_MV.U1q   = KZ04.Uq;
  AT_MV.U1d   = KZ04.Ud;
  AT_LV.U1q   = KZ04.Uq;
  AT_LV.U1d   = KZ04.Ud;
//--------------------------------------(Node 05)-------------------
  AT_MV.U2q = KZ05.Uq;
  AT_MV.U2d = KZ05.Ud;

  L6.U1q = KZ05.Uq;
  L6.U1d = KZ05.Ud;
//--------------------------------------(Node 06)-------------------
  AT_LV.U2q = KZ06.Uq;
  AT_LV.U2d = KZ06.Ud;

  N1.Uqn   = KZ06.Uq;
  N1.Udn   = KZ06.Ud;
//--------------------------------------(Node 07)-------------------
  L1.U2q = KZ07.Uq;
  L1.U2d = KZ07.Ud;
  L6.U2q = KZ07.Uq;
  L6.U2d = KZ07.Ud;

  L3.U1q = KZ07.Uq;
  L3.U1d = KZ07.Ud;
  L4.U1q = KZ07.Uq;
  L4.U1d = KZ07.Ud;

  T3.U1q = KZ07.Uq;
  T3.U1d = KZ07.Ud;
//--------------------------------------(Node 08)-------------------
  T_G1.U2q = KZ08.Uq;
  T_G1.U2d = KZ08.Ud;

  L1.U1q = KZ08.Uq;
  L1.U1d = KZ08.Ud;
//--------------------------------------(Node 09)-------------------
  G1.G.UqG = KZ09.Uq;
  G1.G.UdG = KZ09.Ud;

  T_G1.U1q = KZ09.Uq;
  T_G1.U1d = KZ09.Ud;
//--------------------------------------(Node 10)-------------------
  L4.U2q = KZ10.Uq;
  L4.U2d = KZ10.Ud;
  L5.U2q = KZ10.Uq;
  L5.U2d = KZ10.Ud;

  T1.U1q = KZ10.Uq;
  T1.U1d = KZ10.Ud;
//--------------------------------------(Node 11)-------------------
  T1.U2q = KZ11.Uq;
  T1.U2d = KZ11.Ud;

  N2.Uqn = KZ11.Uq;
  N2.Udn = KZ11.Ud;
//--------------------------------------(Node 12)-------------------
  L3.U2q = KZ12.Uq;
  L3.U2d = KZ12.Ud;

  L5.U1q = KZ12.Uq;
  L5.U1d = KZ12.Ud;

  T2_MV.U1q = KZ12.Uq;
  T2_MV.U1d = KZ12.Ud;
//--------------------------------------(Node 13)-------------------
  T2_HV.U2q = KZ13.Uq;
  T2_HV.U2d = KZ13.Ud;
  T2_MV.U2q = KZ13.Uq;
  T2_MV.U2d = KZ13.Ud;

  T2_LV.U1q = KZ13.Uq;
  T2_LV.U1d = KZ13.Ud;
//--------------------------------------(Node 14)-------------------
  G3.G.UqG = KZ14.Uq;
  G3.G.UdG = KZ14.Ud;

  T2_HV.U1q = KZ14.Uq;
  T2_HV.U1d = KZ14.Ud;
//--------------------------------------(Node 15)-------------------
  T2_LV.U2q = KZ15.Uq;
  T2_LV.U2d = KZ15.Ud;

  N3.Uqn   = KZ15.Uq;
  N3.Udn   = KZ15.Ud;
//--------------------------------------(Node 16)-------------------
  T3.U2q = KZ16.Uq;
  T3.U2d = KZ16.Ud;

  N4.Uqn = KZ16.Uq;
  N4.Udn = KZ16.Ud;
//------------------------------------------------------------------

  L1.Ssys = G1.G.s;
  L2.Ssys = G1.G.s;
  L3.Ssys = G1.G.s;
  L4.Ssys = G1.G.s;
  L5.Ssys = G1.G.s;
  L6.Ssys = G1.G.s;

  N1.Ssys = G1.G.s;
  N2.Ssys = G1.G.s;
  N3.Ssys = G1.G.s;
  N4.Ssys = G1.G.s;
//------------------------------------------------------------------
  G1.G.Tregmax = TregMax;
  G2.G.Tregmax = TregMax;
  G3.G.Tregmax = TregMax;
//-------------------------------------
end Kursovik;
