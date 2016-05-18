clc;
clear all;
%warning off;
rate = 0.0001;
k = 10;
load ex8_movies.mat
A = Y;
[m,n] = size(A);


X = randn(m,k);
Theta = randn(n,k);
v = Theta;
x = X;
lambda = 10;

% r = colperm(A);
% Y = Y(:,r);
% R = R(:,r);


 %spy(R)
% Y = Y';
% R = R';
% r = colperm(Y);
% Y = Y(:,r);
% R = R(:,r);
% Y = Y';
% R = R';

%%
% 
tic
for p = 1:20
    %[J, X] = cofiCostFunc1(X, Theta, Y, R, lambda);
    %J
    [J, X, Theta] = cofiCostFunc2(X, Theta, Y, R, lambda);
    J
end
J_old = J
toc

% x_mat = zeros(m,k);
% t_mat = zeros(n,k);
% ep = 1e-8;
% tic
% for p = 1:100
%     %[J, X] = cofiCostFunc1(X, Theta, Y, R, lambda);
%     %J
%     [J, X_grad, Theta_grad] = cofiCostFunc4(X, Theta, Y, R, lambda);
%     %x_mat = x_mat + X_grad.*X_grad;
%      X = X + rate *X_grad;
%     % t_mat = t_mat + Theta_grad.*Theta_grad;
%     Theta = Theta + rate *Theta_grad;
%     J
% end
% J_old = J
% toc
% idx=find(R(1,:)==1);
%%

% Theta = randn(n,k);
% X = randn(m,k);
% gama1 = 0;
% gama2 = 0.999;
% gama = 0.9;
% iter = 100000;
% x_mat = zeros(m,k);
% x_matrix = zeros(m,k);
% t_mat = zeros(n,k);
% x_mat2 = zeros(m,k);
% t_mat2 = zeros(n,k);
% t_matrix = zeros(n,k);
% ep = 1e-8;
% J = 1;
% tic
% for p = 1:iter
%     [J, X_grad, Theta_grad] = cofiCostFunc(X-x*gama, Theta-v*gama, Y, R, lambda);
% %     v = v*gama + rate*Theta_grad;
% %     x = x*gama + rate*X_grad;
% %     X = X - x;
% %     Theta = Theta - v;
% 
%    %has to have a high initial learning rate = 1
%     x_mat = x_mat + X_grad.*X_grad;
%     X = X - rate ./ sqrt(ep + x_mat).*X_grad;
%     t_mat = t_mat + Theta_grad.*Theta_grad;
%     Theta = Theta - rate ./ sqrt(ep + t_mat).*Theta_grad;
% 
% 
% 
% %     x_mat = x_mat + X_grad.*X_grad;
% %     t_mat = t_mat + Theta_grad.*Theta_grad;
% %  
% %     v = rate ./ sqrt(ep + t_mat).*Theta_grad;
% %     x = rate ./ sqrt(ep + x_mat).*X_grad;
% %     
% %     X = X - x;
% %     Theta = Theta - v;
% 
% 
% 
% %adadelta
% %     x_mat = x_mat + X_grad.*X_grad;
% %     t_mat = t_mat + Theta_grad.*Theta_grad;
% %     
% %     x_matrix = x_matrix + X.*X;
% %     t_matrix = t_matrix + Theta.*Theta;
% %     
% %     %X = X - x;
% %     %Theta = Theta - v;
% %     X = X - sqrt(x_matrix+ep)./sqrt(x_mat+ep).*X_grad; %sqrt(x_matrix+ep)
% %     Theta = Theta - sqrt(t_matrix+ep)./sqrt(t_mat+ep).*Theta_grad;%sqrt(t_matrix+ep)
% 
% %==============================================================%
%     
% 
%     x_mat = gama1*x_mat + (1-gama1)*X_grad; %
%     x_mat2 = gama2*x_mat2 + (1-gama2)*X_grad.*X_grad;
%     m1 = x_mat./(1-gama1^p);
%     m2 = x_mat2 ./(1-gama2^p);
%     X = X - rate ./ (ep + sqrt(m2)).*m1;
%     
%     
%     t_mat = gama1*t_mat + (1-gama1)*Theta_grad; % 
%     t_mat2 = gama2*t_mat2 + (1-gama2)*Theta_grad.*Theta_grad;
%     m1 = t_mat./(1-gama1^p);
%     m2 = t_mat2 ./(1-gama2^p);
%     Theta = Theta - rate ./ (ep + sqrt(m2)).*m1;
% 
% 
% %     X = 
%     %Theta = Theta - rate*Theta_grad;
%      %J
%      J 
%      %if J < J_old
%      %    break;
%      %end
% end
% toc
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % for i = 1:3
% %    [J, X_grad, Theta_grad] = cofiCostFunc3(X, Theta, Y, R, lambda);
% % end
% 
% J
B = X*Theta';
