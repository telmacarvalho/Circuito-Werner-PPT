clear all;
clc;

load('W_parcial_SzOnly_erro.mat');
load('W_Peso.mat');
k=1;
verdadeiro = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
iguais = 0;
for i=1:1001
    contador=0;
    k=i+1;
    if(i==1002)
        disp('Comparação finalizada.')
    else
        for j=k:1001
           resposta= W_parcial_SzOnly_erro(i,:) == W_parcial_SzOnly_erro(j,:);
           contador = contador+1;
           if(resposta==verdadeiro)
               disp('Existem matrizes iguais.')
               n= iguais+1;
               iguais=n;
           end
        end
        if(iguais==0)
            P_Sz(i,1)=1;
        else
            P_Sz(i,1)=0;
        end
        controle(i,1)=contador;
    end    
end

save('P_Sz.mat', 'P_Sz');
