%��������� �������
I_Ug=1; %���������� �� ����� ����������
I_Uc=1; %���������� �� ���
I_Pg=0.8; %�������� �������� �������
Xl=0.1; %������� ������������� � ��
Xt=0.1; %������������� ��
X=Xl-Xt;



%��������� ����������
P_MBt=10; %�������� �������� ����������(���)
Ust_kB=10.5; %����������� ���������� �������(��)
cosFI=0.8; %����������� ����������� �������� ����������
Xd=0.7; %���������� ����������� �������������
Xq=0.49; %���������� ����������� �������������
Xs=0.1; %����������� ������������� ���������
Tj=4; %������������ ����������� ��������
X1d=0.34; %���������� ���������� ����������� �������������
X2d=0.24; %��������������� ���������� ����������� �������������
X2q=0.24; %��������������� ���������� ����������� �������������
Ra=3.236e-3; %�������� ������������� ������� �������
Rr=8.75e-4; %�������� ������������� ������� �����������
Rrd=8.75e-3; %�������� ������������� ����������� ������� � ���������� ���
Rrq=8.75e-3; %�������� ������������� ����������� ������� � ���������� ���
Sigma=0.0475; %����������� �������� ���������� �������
TauC=0.8; %���������� ������� �����������
Tregmax=500; %������ ������������ �������



%��������� ���1
Xl1=0.1; %����������� ������������� ���1
Rl1=0; %�������� ������������� ���1
L1=Xl/100*pi;
Rl=0.001;


%���������
Wc=100*pi;
Mtmax=1.1;

%������
Xad=Xd - Xs;
Xaq=Xq - Xs;
Xsf=1/(1/(X1d - Xs) - 1/Xad);
Xsrd=1/(1/(X2d - Xs) - 1/Xad - 1/Xsf);
Xsrq=1/(1/(X2q - Xs) - 1/Xaq);
X1=(X2d-Xs)/Xsf;
X2=(X2d-Xs)/Xsrd;
X3=(X2q-Xs)/Xsrq;
Xrd=Xsrd + Xad;
Xrq=Xsrq + Xaq;
Xf = Xsf + Xad;
L1=Xl1/Wc;


%������ ������� �������
Base_Sst_BA = P_MBt * 1e6/ cosFI;
Base_Ust_B = Ust_kB * 1e3;
Base_Ist_A = Ust_kB * 1e3;
Base_Srot_BA = Base_Sst_BA;
Base_Irot_A = 2360 * Xad;
Base_Urot_B = Base_Srot_BA * Xad / Base_Irot_A;

%������ ����������� ��������
N_Pg=cosFI;
N_Qg=sqrt(1-(cosFI^2));
N_EQ=sqrt((1 + (Ra*N_Pg + Xq*N_Qg)/1)^2 + ((Xq*N_Pg - Ra*N_Qg)/1)^2);
N_Dg=atan(((Xq*N_Pg - Ra*N_Qg)/1)/(1 + (Ra*N_Pg + Xq*N_Qg)/1));
N_uq=1*cos(N_Dg);
N_ud=-1*sin(N_Dg);
N_id=((N_uq - N_EQ)*Xq - Ra*N_ud)/(Ra^2 + Xq^2);
N_iq=(-N_ud - Ra*N_id)/Xq;
N_Eiq=N_uq - N_id*Xs + Ra*N_iq;
N_ir=-N_id + N_Eiq/Xad;
N_ur=N_ir*Rr;
N_Eq=N_ir*Xad;

%������ ��������� �������
I_Dl=asin((I_Pg*Xl)/(I_Ug*I_Uc));
I_Qg=((I_Ug^2)-I_Ug*I_Uc*cos(I_Dl))/Xl;
I_EQ=sqrt((I_Ug + (Ra*I_Pg + Xq*I_Qg)/I_Ug)^2 + ((Xq*I_Pg - Ra*I_Qg)/I_Ug)^2);
I_Dg=atan(((Xq*I_Pg - Ra*I_Qg)/I_Ug)/(I_Ug + (Ra*I_Pg + Xq*I_Qg)/I_Ug));
I_uq=I_Ug*cos(I_Dg);
I_ud=-I_Ug*sin(I_Dg);
I_id=((I_uq - I_EQ)*Xq - Ra*I_ud)/(Ra^2 + Xq^2);
I_iq=(-I_ud - Ra*I_id)/Xq;
I_Eiq=I_uq - I_id*Xs + Ra*I_iq;
I_ir=-I_id + I_Eiq/Xad;
I_ur=I_ir*Rr;
I_Eq=I_ir*Xad;
I_Mt=I_EQ*I_iq;
I_Mu=I_Mt;
I_Yr=I_Eiq + Xsf*I_ir;
I_Yrq=-I_ud - I_iq*Xs - I_id*Ra;
I_Yrd=I_Eiq;
I_Yd=Xd*I_id + Xad*I_ir;
I_Yq=Xq*I_iq;
I_Dsum = I_Dg + I_Dl;
I_ur_B = I_ur*Base_Urot_B;
I_ir_A = I_ir*Base_Irot_A;

%��������� �������
I_Ug=1; %���������� �� ����� ����������
I_Uc=1; %���������� �� ���
I_Pd=-0.8; %�������� �������� �������
Xl=0.1; %������� ������������� � ��
Xt=0.1; %������������� ��
X=Xl-Xt;



%��������� ���������
P_MBt=10; %�������� �������� ���������(���)
Ust_kB=10.5; %����������� ���������� �������(��)
cosFI=0.8; %����������� ����������� �������� ���������
Xd=0.7; %���������� ����������� �������������
Xq=0.49; %���������� ����������� �������������
Xs=0.1; %����������� ������������� ���������
Tj=4; %������������ ����������� ��������
X1d=0.34; %���������� ���������� ����������� �������������
X2d=0.24; %��������������� ���������� ����������� �������������
X2q=0.24; %��������������� ���������� ����������� �������������
Ra=3.236e-3; %�������� ������������� ������� �������
Rr=8.75e-4; %�������� ������������� ������� �����������
Rrd=8.75e-3; %�������� ������������� ����������� ������� � ���������� ���
Rrq=8.75e-3; %�������� ������������� ����������� ������� � ���������� ���
Sigma=0.0475; %����������� �������� ���������� �������
TauC=0.8; %���������� ������� �����������
Tregmax=500; %������ ������������ �������

%��������� ���1
Xl1=0.1; %����������� ������������� ���1
Rl1=0; %�������� ������������� ���1
L1=Xl/100*pi;
Rl=0.001;


%���������
Wc=100*pi;
Mtmax=1.1;

%������
Xad=Xd - Xs;
Xaq=Xq - Xs;
Xsf=1/(1/(X1d - Xs) - 1/Xad);
Xsrd=1/(1/(X2d - Xs) - 1/Xad - 1/Xsf);
Xsrq=1/(1/(X2q - Xs) - 1/Xaq);
X1=(X2d-Xs)/Xsf;
X2=(X2d-Xs)/Xsrd;
X3=(X2q-Xs)/Xsrq;
Xrd=Xsrd + Xad;
Xrq=Xsrq + Xaq;
Xf = Xsf + Xad;
L1=Xl1/Wc;


%������ ������� �������
Base_Sst_BA = P_MBt * 1e6/ cosFI;
Base_Ust_B = Ust_kB * 1e3;
Base_Ist_A = Ust_kB * 1e3;
Base_Srot_BA = Base_Sst_BA;
Base_Irot_A = 2360 * Xad;
Base_Urot_B = Base_Srot_BA * Xad / Base_Irot_A;

%������ ����������� ��������
N_Pg=cosFI;
N_Qg=sqrt(1-(cosFI^2));
N_EQ=sqrt((1 + (Ra*N_Pg + Xq*N_Qg)/1)^2 + ((Xq*N_Pg - Ra*N_Qg)/1)^2);
N_Dg=atan(((Xq*N_Pg - Ra*N_Qg)/1)/(1 + (Ra*N_Pg + Xq*N_Qg)/1));
N_uq=1*cos(N_Dg);
N_ud=-1*sin(N_Dg);
N_id=((N_uq - N_EQ)*Xq - Ra*N_ud)/(Ra^2 + Xq^2);
N_iq=(-N_ud - Ra*N_id)/Xq;
N_Eiq=N_uq - N_id*Xs + Ra*N_iq;
N_ir=-N_id + N_Eiq/Xad;
N_ur=N_ir*Rr;
N_Eq=N_ir*Xad;

%������ ��������� �������
I_Dl_2=asin((I_Pd*Xl)/(I_Ug*I_Uc));
I_Qg_2=((I_Ug^2)-I_Ug*I_Uc*cos(I_Dl_2))/Xl;
I_EQ_2=sqrt((I_Ug + (Ra*I_Pd + Xq*I_Qg_2)/I_Ug)^2 + ((Xq*I_Pd - Ra*I_Qg)/I_Ug)^2);
I_Dg_2=atan(((Xq*I_Pd - Ra*I_Qg_2)/I_Ug)/(I_Ug + (Ra*I_Pd + Xq*I_Qg)/I_Ug));
I_uq_2=I_Ug*cos(I_Dg_2);
I_ud_2=-I_Ug*sin(I_Dg_2);
I_id_2=((I_uq_2 - I_EQ_2)*Xq - Ra*I_ud_2)/(Ra^2 + Xq^2);
I_iq_2=(-I_ud_2 - Ra*I_id_2)/Xq;
I_Eiq_2=I_uq_2 - I_id_2*Xs + Ra*I_iq_2;
I_ir_2=-I_id_2 + I_Eiq_2/Xad;
I_ur_2=I_ir_2*Rr;
I_Eq=I_ir_2*Xad;
I_Mt_2=I_EQ_2*I_iq_2;
I_Mu_2=I_Mt_2;
I_Yr_2=I_Eiq_2 + Xsf*I_ir_2;
I_Yrq_2=-I_ud - I_iq_2*Xs - I_id_2*Ra;
I_Yrd_2=I_Eiq_2;
I_Yd_2=Xd*I_id_2 + Xad*I_ir_2;
I_Yq_2=Xq*I_iq_2;
I_Dsum_2 = I_Dg_2 + I_Dl_2;
I_ur_B_2 = I_ur_2*Base_Urot_B;
I_ir_A_2 = I_ir_2*Base_Irot_A;