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

p=0.5;

rho=p*rho_bell1+(1-p)*rho_bell2;

Traco = rho(1,1)+rho(2,2)+rho(3,3)+rho(4,4);    
assert(Traco < (1+exp(-10)) && Traco > (1-exp(-10)))