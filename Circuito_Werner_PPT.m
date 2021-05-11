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

Parametro = -0.01;
for n=1:101
    % Estabelecendo os pesos
    P1 = Parametro+0.01;
    P2 = (1-P1)/4;
    Entrada = {sqrt(P1)*Bell1 (sqrt(P2))*E0 (sqrt(P2))*E01 (sqrt(P2))*E10...
        (sqrt(P2))*E11};
    W_P_101(n,1) = P1;
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
    
    % Tomografia unitária
   Rho_parcial_SzSy = Tomografia_parcial_SzSy(Rhos);

   %Armazenamento da matriz densidade parcial
   Werner_parcial = reshape(Rho_parcial_SzSy, 1, 16);
   Werner_parcial1 = real(Werner_parcial);
   W_parcial_101(n, 1:16) = (Werner_parcial1);
        
    % Armazenamento da matriz densidade do estado de Werner  
    % Conversão de dados
    Werner = reshape(Rhos, 1, 16);
    Werner1 = real(Werner);
    W_estados_101(n,1:16) = (Werner1);
        
    
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
        x(n) = W_P_101(n);
        y(n) = 1;
    else
        x(n) = -1;
        y(n) = -1;
    end
    if (Resultado(n) == 0)
        z(n) = W_P_101(n);
        k(n) = 0;
    else
        z(n) = -1;
        k(n) = -1;
    end
end

% Gráfico em barra 
%x = [P];
%y = [Resultado];
%figure(1), clf%
%bar (x,y)
%xlabel('Peso (P)')
%ylabel('Estados separáveis')
%title('Estados separáveis pelo critério PPT por peso')

% Exportação de dados
Rotulos = categorical(Rotulos);
W_PPT_101 = dummyvar(Rotulos);
Dados = horzcat(W_estados_101, W_PPT_101);
%csvwrite('Dados_1001.csv', Dados);
save('W_estados_101.mat', 'W_estados_101');
save('W_PPT_101.mat','W_PPT_101');
save('W_P_101.mat', 'W_P_101');
save('W_parcial_101.mat', 'W_parcial_101');

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
title('Classificação por peso utilizando o critério PPT')
% Fixando um sombreado para P>1/3
cm = [0 0 0 ;  0.9 0.9 0.9;  1 1 1];
patch([(1/3) (1/3) 1.01 1.01 (1/3)]', [-0.95 1.95 1.95 -0.95 -0.95]', cm(2,:), 'EdgeColor','none', 'DisplayName', 'Área de emaranhamento')
axis([-0.02 1.02 -1 2])
set(gca,'children',flipud(get(gca,'children')))
