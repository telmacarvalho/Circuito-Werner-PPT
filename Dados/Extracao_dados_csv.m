% Limpeza da command window e de dados
clc;
clear all;

% Importação dos dados
load('B_completo.mat')
load('B_PPT.mat')
load('B_parcial_SzSy_erro.mat')
load('B_parcial_SzSx_erro.mat')
load('B_parcial_SzOnly_erro.mat')

% Mudando o formato dos dados
save('B_completo.csv', 'B_completo');
save('B_PPT.csv','B_PPT');
save('B_parcial_SzSy_erro.csv', 'B_parcial_SzSy_erro');
save('B_parcial_SzSx_erro.csv', 'B_parcial_SzSx_erro');
save('B_parcial_SzOnly_erro.csv', 'B_parcial_SzOnly_erro');