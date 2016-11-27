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

function Lin_Law
  input Real x1, x2, y1, y2, x;
  output Real y;
  Real a, b;
algorithm
  a = (y1-y2)/(x1-x2);
  b = y1 - a*x1;
  y = a*x + b;
end Lin_Law;

model ExcitationRegulator "Simple model of ARV-SDP1"
  parameter Real K0u=-10;
  parameter Real K1u=0;
  parameter Real K0w=0;
  parameter Real K1w=0;
  parameter Real K1if=0;
  parameter Real DL0=0;

  parameter Real Tokp=0.05;

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
  parameter Real Gn=0.05 "Активная составляющая нагрузки";
  parameter Real Bn=-0.025 "Реактивная составляющая нагрузки";
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
  Gnn = Gn; //             R0/(R0^2 + (X0*(1+Ssys))^2);
  Bnn = Bn; // -(X0*(1+Ssys))/(R0^2 + (X0*(1+Ssys))^2);
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

model HVLine "The model of transmission line"
//---------------------------------------------
  parameter Real Rline=0.01;
  parameter Real Xline=0.1;
  parameter Real Bline=0;
//---------------------------------------------
  parameter Real TLineOff=1000;
  parameter Real dTLineOff=10;
  parameter Real Koff=100000;
//---------------------------------------------
  Real U1d, U1q, U2d, U2q, I12d, I12q, I12m, U1m, U2m;
  Real I1d, I1q, I2d, I2q, Idc1, Iqc1, Idc2, Iqc2, DU1, DU2, Ssys, Bnn;
  Real P1, Q1, P2, Q2, RL, XL, BL;
equation
// the algorithm of triping transmission line
  RL = if time>=TLineOff and time<TLineOff+dTLineOff then Koff*Rline
                                                     else Rline;
  XL = if time>=TLineOff and time<TLineOff+dTLineOff then Koff*Xline
                                                     else Xline;
  BL = if time>=TLineOff and time<TLineOff+dTLineOff then Bline/Koff
                                                     else Bline;
// series impedance between nodes
  U1q = U2q + RL*I12q - XL*I12d;
  U1d = U2d + RL*I12d + XL*I12q;
// line capacitances
  Bnn = BL/(1+Ssys);
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

model SyncronousMachine "Complete model of synchronous machine" 
//------------------------------------------------------------------------------------
  parameter Real TgenOff=1000;
  parameter Real dTgenOff=10;
//------------------------------------------------------------------------------------
  parameter Real Pg=-0.35;
  parameter Real Qg=-0.2;
  parameter Real Ut=1;
//------------------------------------------------------------------------------------
  parameter Real Xd=3.223 "Продольное индуктивное сопротивление";
  parameter Real Xq=2.089 "Поперечное индуктивное сопротивление";
  parameter Real Xs=0.257 "Индуктивное сопротивление рассеяния";
  parameter Real X1d=0.829 "Переходное продольное индуктивное сопротивление";
  parameter Real X2d=0.485 "Сверхпереходное продольное индуктивное сопротивление";
  parameter Real X2q=0.457 "Сверхпереходное поперечное индуктивное сопротивление";
  parameter Real Rf=0.00178 "Активное сопротивление обмотки возбуждения";
  parameter Real R1d=0.0127 "Активное сопротивление ДК в продольной оси";
  parameter Real R1q=0.00871 "Активное сопротивление ДК в поперечной оси";
  parameter Real Tj=2.45 "Механическая инерционная постоянная";
  parameter Real Xt=0 "Индуктивное сопротивление трансформатора";
  parameter Real Sigma=0.0475 "Коэффициент статизма регулятора турбины";
  parameter Real TauC=99999 "Постоянная времени сервомотора";

  parameter Real Mt_max=0 "Верхнее ограничение момента турбины";
  parameter Real Mt_min=-0.5   "Нижнее ограничение момента турбины";
//------------------------------------------------------------------------------------
  Real dUtr, dWu, dIf, dEr, dWf;
  constant Real PI=4*arctan(1);
  constant Real Wc=100*PI;
// Parameters of power system elements
  constant Real Ra=0;
  constant Real Rt=0;
  //------------------------------------------------------------------
  Real EQ, Dg, Uq, Ud, Id, Iq, Mt, Eiq, Ir, Uf;
  Real Eq, Er0, PsiR, PsiRD, PsiRQ, Er_max, Er_min, Uf_full;
// Uc, Dl, Ug, Dt, Pt, Qt, , X2ds, X2qs, Udgen, Uqgen, Ugen0,
  Real Xad, Xaq, X2dt, X2qt, Rs, Xsf, Xs1d, Xs1q, X1, X2, X3;
  Real Ugen;
// Ut, Pg, Qg, 
  Real Mu0, Ro, Mu, Mt_pp;
  Real Pgen, Qgen;
  Real UdG, UqG, IdG, IqG, DeltaIJ, Ssys;
  //------------------------------------------------------------------
// Integrated Variables
  Real Yr, Yrd, Yrq, s;
//, Dsum
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
    Dg=arctan2((((Xq+Xt)*Pg-(Ra+Rt)*Qg)/Ut),(Ut+((Ra+Rt)*Pg+(Xq+Xt)*Qg)/Ut));
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
//    reinit(Dsum, Dg);
    reinit(s, 0);
    reinit(Yr, PsiR);
    reinit(Yrd, PsiRD);
    reinit(Yrq, PsiRQ);
    reinit(Mu, Mu0);
  end when;
//  der(Dsum) = Wc*s;

  Uf_full = if      (Uf*Xad/Rf+dEr) > Er_max then Er_max*Rf/Xad
            else if (Uf*Xad/Rf+dEr) < Er_min then Er_min*Rf/Xad
                                             else Uf + Rf*dEr/Xad;

// or ((time>Tkz_R1) and (time<Tkz_R1+dTforce)))
// if (((time>Tkz_G) and (time<Tkz_G+dTforce)) then Er_max*Rf/Xad 

  der(Yr) = Wc*(Uf_full - Rf*ir);
  der(Yrd) = -Wc*R1d*ird;
  der(Yrq) = -Wc*R1q*irq;
  der(s) =(Mt_pp - Me)/Tj;
  Ro = -Mu + Mu0 - s/Sigma;
  der(Mu) = Ro/TauC;

  Mt_pp = if Mu>=Mt_min and Mu<=Mt_max then Mu
     else if Mu>Mt_max                 then Mt_max
                                       else Mt_min;

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

model vankor11 "The simplified model of autonomous power system"
  Generator_with_ARV G1;

  HVLine L0102;

  Generator_with_ARV SM2;

  ShortCircuitShunt KZ1;
  ShortCircuitShunt KZ2;


  Real delta_G1_SM2;

equation
//-------------------------The equations of current balance----------------------------------
//--------------------------------------------(Node 1)-------------------------------------------------
  G1.G.IdG - L0102.I1d - KZ1.id = 0;
  G1.G.IqG - L0102.I1q - KZ1.iq = 0;
//--------------------------------------------(Node 2)-------------------------------------------------
  SM2.G.IdG + L0102.I2d - KZ2.id = 0;
  SM2.G.IqG + L0102.I2q - KZ2.iq = 0;

//  SM2.G.IdG + IM2.Id + L0102.I2d - N2.Idn - R2.Idn - KZ2.id = 0;
//  SM2.G.IqG + IM2.Iq + L0102.I2q - N2.Iqn - R2.Iqn - KZ2.iq = 0;
//-----------------------------------------------------------------------------------------------------


//-----------Angles of Machines---------------------------
  der(delta_G1_SM2) = G1.G.Wc*(G1.G.s - SM2.G.s);
  
  G1.G.DeltaIJ  = 0;
  SM2.G.DeltaIJ = delta_G1_SM2;

  G1.G.Ssys  = G1.G.s;
  SM2.G.Ssys = G1.G.s;


//---------------------------Line connections----------------------
//--------------------------------------(Node 1)-------------------
  G1.G.UqG = KZ1.Uq;
  G1.G.UdG = KZ1.Ud;

  L0102.U1q   = KZ1.Uq;
  L0102.U1d   = KZ1.Ud;
//--------------------------------------(Node 2)-------------------
  SM2.G.UqG = KZ2.Uq;
  SM2.G.UdG = KZ2.Ud;

//  IM2.Uq = KZ2.Uq;
//  IM2.Ud = KZ2.Ud;

  L0102.U2q = KZ2.Uq;
  L0102.U2d = KZ2.Ud;

//  N2.Uqn = KZ2.Uq;
//  N2.Udn = KZ2.Ud;

//  R2.Uqn = KZ2.Uq;
//  R2.Udn = KZ2.Ud;
//-----------------------------------------------------------------


  L0102.Ssys = G1.G.s;

//  N2.Ssys = G1.G.s;
//  R2.Ssys = G1.G.s;

//  R2.ILine = L0102.I12m;
//-------------------------------------
end vankor11;
