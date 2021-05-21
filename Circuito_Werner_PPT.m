% Limpeza da command window e de dados
clc;
clear all;


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
E0 = [0; 0; 0; 0];

% Entrada de estados
E00 = HH;
E01 = HV;
E10 = VH;
E11 = VV;
Bell1 = 1/sqrt(2)*(kron(H,H) + kron(V,V));
Bell2 = 1/sqrt(2)*(kron(H,V) + kron(V,H));
Bell3 = 1/sqrt(2)*(kron(H,H) - kron(V,V));
Bell4 = 1/sqrt(2)*(kron(H,V) - kron(V,H));

% Estabelecendo quantidade de estados gerados 
Quant_dados = 101;

if (Quant_dados == 101)
    N = 101;
    Parametro = -0.01;
elseif (Quant_dados == 1001)
    N = 1001;
    Parametro = -0.001;
elseif (Quant_dados == 10001)
    N = 10001;
    Parametro = -0.0001;
else
    N = 11;
    Parametro = -0.1;
end

for n=1:N
    % Estabelecendo os pesos
    if (Quant_dados == 101)
         P1 = Parametro+0.01;
    elseif (Quant_dados == 1001)
        P1 = Parametro+0.001;
    elseif (Quant_dados == 10001)
        P1 = Parametro+0.0001;
    else
        P1 = Parametro+0.1;
    end   
    P2 = (1-P1)/4;
    Entrada = {sqrt(P1)*Bell1 (sqrt(P2))*E0 (sqrt(P2))*E01 (sqrt(P2))*E10...
        (sqrt(P2))*E11};
    W_Peso(n,1) = P1;
    Parametro = P1;
    Rhos = 0;


    % Tomografia dos 5 estados    
    for m=1:5
        Rho = Tomografia(Entrada{m});        
        % Soma das matrizes densidade
        Soma = Rhos+Rho;
        Rhos = Soma;
    end

    % Teste do traço unitário
    Traco = Rhos(1,1)+Rhos(2,2)+Rhos(3,3)+Rhos(4,4);    
    assert(Traco < (1+exp(-10)) && Traco > (-1-exp(-10)))

    % Tomografia parcial
    Rho_parcial_SzSy = Tomografia_parcial_SzSy_erro(Rhos); % Não foi medido em Sx
    Rho_parcial_SzSx = Tomografia_parcial_SzSx_erro(Rhos); % Não foi medido em Sy
    Rho_parcial_SzSxParcial = Tomografia_parcial_SzSxParcial_erro(Rhos); % Não foi medido em...
    % Sy e mediu parcialmente Sx

    %Armazenamento da matriz densidade parcial
    Werner_parcial_SzSy = reshape(Rho_parcial_SzSy, 1, 16);
    Werner_parcial1 = real(Werner_parcial_SzSy);
    W_parcial_SzSy_erro(n, 1:16) = (Werner_parcial1);


    Werner_parcial_SzSx = reshape(Rho_parcial_SzSx, 1, 16);
    Werner_parcial2 = real(Werner_parcial_SzSx);
    W_parcial_SzSx_erro(n, 1:16) = (Werner_parcial2);
    
    Werner_parcial_SzSxParcial = reshape(Rho_parcial_SzSxParcial, 1, 16);
    Werner_parcial3 = real(Werner_parcial_SzSxParcial);
    W_parcial_SzSxParcial_erro(n, 1:16) = (Werner_parcial3);

    % Armazenamento da matriz densidade do estado de Werner  
    % Conversão de dados
    Werner = reshape(Rhos, 1, 16);
    Werner1 = real(Werner);
    W_completo(n,1:16) = (Werner1);


    % Cálculo PPT
    % Transposição parcial em relação ao sistema A
    a1 = [Rhos(1,3:4); Rhos(2,3:4)];
    a2 = [Rhos(3,1:2); Rhos(4,1:2)];
    A1 = a1';
    A2 = a2';
    B1 = [Rhos(1,1:2); Rhos(2,1:2)];
    B2 = [Rhos(3,3:4); Rhos(4,3:4)];
    RhosPPT = [B1 A1; A2 B2];

    % Calculando os autovalores
    Autovalores = eig(RhosPPT);

    % Definindo se o estado é emaranhado: emaranhado = 0 e separável = 1
    if (Autovalores(1)<(-1*exp(-10)))
        AutNeg=0;
        Rotulo = 'Emaranhado';
    elseif (Autovalores(2)<(-1*exp(-10)))
        AutNeg=0;
        Rotulo = 'Emaranhado';
    elseif (Autovalores(3)<(-1*exp(-10)))
        AutNeg=0;
        Rotulo = 'Emaranhado';
    elseif (Autovalores(4)<(-1*exp(-10)))
        AutNeg=0;
        Rotulo = 'Emaranhado';
    else
        AutNeg=1;
        Rotulo = 'Separável';
    end

    %  Teste do resultado: estado de Werner é emaranho se P1 > 1/3
    if(P1 > 1/3)
        assert(AutNeg == 0)
    else
        assert(AutNeg == 1)
    end

    % Armazenamento do resultado PPT
    Resultado(n,1) = (AutNeg);
    Rotulos{n,1} = (Rotulo);

    % Preparando dados para o gráfico classificatório
    if (Resultado(n) == 1)
        x(n) = W_Peso(n);
        y(n) = 1;
    else
        x(n) = -1;
        y(n) = -1;
    end
    if (Resultado(n) == 0)
        z(n) = W_Peso(n);
        k(n) = 0;
    else
        z(n) = -1;
        k(n) = -1;
    end
end

% Exportação de dados
Rotulos = categorical(Rotulos);
W_PPT = dummyvar(Rotulos);
save('W_completo.mat', 'W_completo');
save('W_PPT.mat','W_PPT');
save('W_Peso.mat', 'W_Peso');
save('W_parcial_SzSy_erro.mat', 'W_parcial_SzSy_erro');
save('W_parcial_SzSx_erro.mat', 'W_parcial_SzSx_erro');
save('W_parcial_SzSxParcial_erro.mat', 'W_parcial_SzSxParcial_erro');


% Gráfico classificatório
x = x(x>=0);
y = y(y>=0);
z = z(z>=0);
k = k(k>=0);
figure
plot(x, y, 'b.', z, k, 'r.', 'MarkerSize', 20)
set(gca,'FontSize',18)
set(gca, 'FontName', 'Times New Roman'); 
%set(gca,'Color','none')
xticks([0:0.1:1])
yticks([0 1])
yticklabels({ })
legend({'Separável','Emaranhado'},'Location','southwest', 'Color','none')
title('Classificação por peso utilizando o critério PPT (dados com erro nos ângulos)')
% Fixando um sombreado para P>1/3
cm = [0 0 0 ;  0.9 0.9 0.9;  1 1 1];
patch([(1/3) (1/3) 1.01 1.01 (1/3)]', [-0.95 1.95 1.95 -0.95 -0.95]', cm(2,:), 'EdgeColor','none', 'DisplayName', 'Área de emaranhamento')
axis([-0.02 1.02 -1 2])
set(gca,'children',flipud(get(gca,'children')))
