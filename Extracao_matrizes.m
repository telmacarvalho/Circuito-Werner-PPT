
%load('W_completo.mat')

for n = 1:1001
    %Werner = reshape(W_completo(n, 1:16), 4, 4);
    Werner = W_completo(n, 1:16);
    Werner1 = reshape(Werner, 4, 4);
    W_original{n,1} = Werner1;
end
save('W_original.mat', 'W_original');
