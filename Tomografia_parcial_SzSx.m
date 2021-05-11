function [Rho] = Tomografia(Entrada)
% Dispositivos ópticos
% Variáveis importantes
x1 = pi/8;
x2 = 0;
y1 = pi/8;
y2 = 0;
y3 = -pi/2;

%HWP
HWP = [[cos(2*x1), sin(2*x1)]; [sin(2*x1), -cos(2*x1)]];

%QWP
QWP = [[cos(2*x2)^2-1i*(sin(2*x2))^2, cos(2*x2)*sin(2*x2)*(1+1i)];...
    [cos(2*x2)*sin(2*x2)*(1+1i), -1i*(cos(2*x2))^2+(sin(2*x2))^2]];

%DP
DP = [[cos(2*y1), sin(2*y1)]; [sin(2*y1), -cos(2*y1)]];

%CL
CL = [[cos(y2)^2+(exp(1i*y3))*(sin(y2))^2, cos(y2)*sin(y2)*(exp(1i*y3)-1)];...
    [(exp(1i*y3)-1)*sin(y2)*cos(y2), exp(1i*y3)*(cos(y2))^2+(sin(y2))^2]];

%Beam Splitter matrix 50/50
BS = 1/sqrt(2)*[1 1; 1 -1];


%Matrizes básicas
%Matriz identidade
I = [1 0; 0 1];
%Matrizes de Pauli
s0 = I;
s1 = [0 1;1 0];
s2 = [0 -1i; 1i 0];
s3 = [1 0; 0 -1];

% Bases de medida
H = [1;0];
V = [0;1];
D = 1/sqrt(2)*[1;1];
A = 1/sqrt(2)*[1;-1];
R = 1/sqrt(2)*[1;-1i];
L = 1/sqrt(2)*[1;1i];


%Variáveis para projeções
HH = kron(H,H);
HV = kron(H,V);
VH = kron(V,H);
VV = kron(V,V);

DD = kron(D,D);
DA = kron(D,A);
AD = kron(A,D);
AA = kron(A,A);

RR = kron(R,R);
RL = kron(R,L);
LR = kron(L,R);
LL = kron(L,L);

HD = kron(H,D);
HA = kron(H,A);
HR = kron(H,R);
HL = kron(H,L);

VD = kron(V,D);
VA = kron(V,A);
VR = kron(V,R);
VL = kron(V,L);

RD = kron(R,D);
RA = kron(R,A);
RV = kron(R,V);
RH = kron(R,H);

LD = kron(L,D);
LA = kron(L,A);
LV = kron(L,V);
LH = kron(L,H);

DH = kron(D,H);
DV = kron(D,V);
DR = kron(D,R);
DL = kron(D,L);

AH = kron(A,H);
AV = kron(A,V);
AR = kron(A,R);
AL = kron(A,L);

% Circuito: probabilidades
PHH = ((HH)'*Entrada*HH)'*((HH)'*Entrada*HH);
PHV = ((HV)'*Entrada*HV)'*((HV)'*Entrada*HV);
PVH = ((VH)'*Entrada*VH)'*((VH)'*Entrada*VH);
PVV = ((VV)'*Entrada*VV)'*((VV)'*Entrada*VV);

PHD = ((HH)'*kron(I,DP)*Entrada*HH)'*((HH)'*kron(I,DP)*Entrada*HH);
PHA = ((HV)'*kron(I,DP)*Entrada*HV)'*((HV)'*kron(I,DP)*Entrada*HV);
PVD = ((VH)'*kron(I,DP)*Entrada*VH)'*((VH)'*kron(I,DP)*Entrada*VH);
PVA = ((VV)'*kron(I,DP)*Entrada*VV)'*((VV)'*kron(I,DP)*Entrada*VV);

PHR = ((HH)'*kron(I,DP*CL)*Entrada*HH)'*((HH)'*kron(I,DP*CL)*Entrada*HH);
PHL = ((HV)'*kron(I,DP*CL)*Entrada*HV)'*((HV)'*kron(I,DP*CL)*Entrada*HV);
PVR = ((VH)'*kron(I,DP*CL)*Entrada*VH)'*((VH)'*kron(I,DP*CL)*Entrada*VH);
PVL = ((VV)'*kron(I,DP*CL)*Entrada*VV)'*((VV)'*kron(I,DP*CL)*Entrada*VV);

PDD = ((HH)'*kron(HWP,DP)*Entrada*HH)'*((HH)'*kron(HWP,DP)*Entrada*HH);
PDA = ((HV)'*kron(HWP,DP)*Entrada*HV)'*((HV)'*kron(HWP,DP)*Entrada*HV);
PAD = ((VH)'*kron(HWP,DP)*Entrada*VH)'*((VH)'*kron(HWP,DP)*Entrada*VH);
PAA = ((VV)'*kron(HWP,DP)*Entrada*VV)'*((VV)'*kron(HWP,DP)*Entrada*VV);

PDH = ((HH)'*kron(HWP,I)*Entrada*HH)'*((HH)'*kron(HWP,I)*Entrada*HH);
PDV = ((HV)'*kron(HWP,I)*Entrada*HV)'*((HV)'*kron(HWP,I)*Entrada*HV);
PAH = ((VH)'*kron(HWP,I)*Entrada*VH)'*((VH)'*kron(HWP,I)*Entrada*VH);
PAV = ((VV)'*kron(HWP,I)*Entrada*VV)'*((VV)'*kron(HWP,I)*Entrada*VV);

PDR = ((HH)'*kron(HWP,DP*CL)*Entrada*HH)'*((HH)'*kron(HWP,DP*CL)*Entrada*HH);
PDL = ((HV)'*kron(HWP,DP*CL)*Entrada*HV)'*((HV)'*kron(HWP,DP*CL)*Entrada*HV);
PAL = ((VH)'*kron(HWP,DP*CL)*Entrada*VH)'*((VH)'*kron(HWP,DP*CL)*Entrada*VH);
PAR = ((VV)'*kron(HWP,DP*CL)*Entrada*VV)'*((VV)'*kron(HWP,DP*CL)*Entrada*VV);

PRR = ((HH)'*kron(HWP*QWP,DP*CL)*Entrada*HH)'*((HH)'*kron(HWP*QWP,DP*CL)*Entrada*HH);
PRL = ((HV)'*kron(HWP*QWP,DP*CL)*Entrada*HV)'*((HV)'*kron(HWP*QWP,DP*CL)*Entrada*HV);
PLR = ((VH)'*kron(HWP*QWP,DP*CL)*Entrada*VH)'*((VH)'*kron(HWP*QWP,DP*CL)*Entrada*VH);
PLL = ((VV)'*kron(HWP*QWP,DP*CL)*Entrada*VV)'*((VV)'*kron(HWP*QWP,DP*CL)*Entrada*VV);

PRH = ((HH)'*kron(HWP*QWP,I)*Entrada*HH)'*((HH)'*kron(HWP*QWP,I)*Entrada*HH);
PRV = ((HV)'*kron(HWP*QWP,I)*Entrada*HV)'*((HV)'*kron(HWP*QWP,I)*Entrada*HV);
PLH = ((VH)'*kron(HWP*QWP,I)*Entrada*VH)'*((VH)'*kron(HWP*QWP,I)*Entrada*VH);
PLV = ((VV)'*kron(HWP*QWP,I)*Entrada*VV)'*((VV)'*kron(HWP*QWP,I)*Entrada*VV);

PRD = ((HH)'*kron(HWP*QWP,DP)*Entrada*HH)'*((HH)'*kron(HWP*QWP,DP)*Entrada*HH);
PRA = ((HV)'*kron(HWP*QWP,DP)*Entrada*HV)'*((HV)'*kron(HWP*QWP,DP)*Entrada*HV);
PLD = ((VH)'*kron(HWP*QWP,DP)*Entrada*VH)'*((VH)'*kron(HWP*QWP,DP)*Entrada*VH);
PLA = ((VV)'*kron(HWP*QWP,DP)*Entrada*VV)'*((VV)'*kron(HWP*QWP,DP)*Entrada*VV);

% Produto tensorial entre as matrizes de Pauli
s00 = kron(s0,s0);
s01 = kron(s0,s1);
s02 = kron(s0,s2);
s03 = kron(s0,s3);
s10 = kron(s1,s0);
s11 = kron(s1,s1);
s12 = kron(s1,s2);
s13 = kron(s1,s3);
s20 = kron(s2,s0);
s21 = kron(s2,s1);
s22 = kron(s2,s2);
s23 = kron(s2,s3);
s30 = kron(s3,s0);
s31 = kron(s3,s1);
s32 = kron(s3,s2);
s33 = kron(s3,s3);

% Parâmetros de Stokes em termos de probabilidades
S00 = PHH+PHV+PVH+PVV;
S01 = PHD-PHA+PVD-PVA;
S02 = PHR-PHL+PVR-PVL;
S03 = PHH-PHV+PVH-PVV;

S10 = PDH+PDV-PAH-PAV;
S11 = PDD-PDA-PAD+PAA;
S12 = PDR-PDL-PAR+PAL;
S13 = PDH-PDV-PAH+PAV;

S20 = PRH+PRV-PLH-PLV;
S21 = PRD-PRA-PLD+PLA;
S22 = PRR-PRL-PLR+PLL;
S23 = PRH-PRV-PLH+PLV;

S30 = PHH+PHV-PVH-PVV;
S31 = PHD-PHA-PVD+PVA;
S32 = PHR-PHL-PVR+PVL;
S33 = PHH-PHV-PVH+PVV;

% Matriz densidade calculada através das probalidades de projeção
Rho = (1/4)*(S00*s00+S01*s01+S02*s02+S03*s03+S10*s10+S11*s11+S12*s12+...
    S13*s13+S20*s20+S21*s21+S22*s22+S23*s23+S30*s30+S31*s31+S32*s32+S33*s33);

% Para conferir o resultado da matriz densidade
rho = Entrada*(Entrada)';
    for j=1:16
        assert(Rho(j) < (rho(j)+exp(-10)) && Rho(j) > (rho(j)-exp(-10)))
    end
end