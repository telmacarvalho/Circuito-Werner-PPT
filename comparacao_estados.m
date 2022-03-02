clear all;
clc;

load('W_parcial_SzOnly_erro.mat');
k=1;
verdadeiro = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];

for i=1:1001
    n=0;
    k=i+1;
    if(i==1001)
        disp('Comparação finalizada.')
    else
        for j=k:1001
           resposta= W_parcial_SzOnly_erro(i,:) == W_parcial_SzOnly_erro(j,:);
           if(resposta==verdadeiro)
               disp('Existem matrizes iguais.')
           else
           n = n+1;
           end
        end
        testes(i,1)=n;
    end    
end