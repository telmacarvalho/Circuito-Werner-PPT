clear all;
clc;

H = [1;0];
V = [0;1];

HH = kron(H,H);
HV = kron(H,V);
VH = kron(V,H);
VV = kron(V,V);

E00 = HH;
E01 = HV;
E10 = VH;
E11 = VV;
Bell1 = 1/sqrt(2)*(kron(H,H) + kron(V,V));
Bell2 = 1/sqrt(2)*(kron(H,V) + kron(V,H));
Bell3 = 1/sqrt(2)*(kron(H,H) - kron(V,V));
Bell4 = 1/sqrt(2)*(kron(H,V) - kron(V,H));

rho_bell1=Bell1*(Bell1)';
rho_bell2=Bell2*(Bell2)';
rho_bell3=Bell3*(Bell3)';
rho_bell4=Bell4*(Bell4)';

p=0.5;

rho1=p*rho_bell1+(1-p)*rho_bell2;
rho2=p*rho_bell3+(1-p)*rho_bell4;

p1=0.5;

Rho= p1*rho1+(1-p1)*rho2;

Traco2 = rho2(1,1)+rho2(2,2)+rho2(3,3)+rho2(4,4);    
assert(Traco2 < (1+exp(-10)) && Traco2 > (1-exp(-10)))

Traco = Rho(1,1)+Rho(2,2)+Rho(3,3)+Rho(4,4);    
assert(Traco < (1+exp(-10)) && Traco > (1-exp(-10)))