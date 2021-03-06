function [J, X, Theta] = cofiCostFunc2(X, Theta, Y, R, lambda)
%COFICOSTFUNC Collaborative filtering cost function
%   [J, grad] = COFICOSTFUNC(params, Y, R, num_users, num_movies, ...
%   num_features, lambda) returns the cost and gradient for the
%   collaborative filtering problem.

            
% You need to return the following values correctly
J = 0;
[m,n] = size(Y);
% ====================== YOUR CODE HERE ======================
% Instructions: Compute the cost function and gradient for collaborative
%               filtering. Concretely, you should first implement the cost
%               function (without regularization) and make sure it is
%               matches our costs. After that, you should implement the 
%               gradient and use the checkCostFunction routine to check
%               that the gradient is correct. Finally, you should implement
%               regularization.
%
% Notes: X - num_movies  x num_features matrix of movie features
%        Theta - num_users  x num_features matrix of user features
%        Y - num_movies x num_users matrix of user ratings of movies
%        R - num_movies x num_users matrix, where R(i, j) = 1 if the 
%            i-th movie was rated by the j-th user
%
% You should set the following variables correctly:
%
%        X_grad - num_movies x num_features matrix, containing the 
%                 partial derivatives w.r.t. to each element of X
%        Theta_grad - num_users x num_features matrix, containing the 
%                     partial derivatives w.r.t. to each element of Theta
%
%J = 0;

J=1/2*sum(sum(R.*((X*Theta'-Y).^2)))+lambda/2*sum(sum(Theta.^2))+lambda/2*sum(sum(X.^2));

% for i=1:m
%   idx=find(R(i,:)==1);
%   theta_tmp=(Theta(idx,:))';
%   A = zeros(size(X,2));
%   A = A + theta_tmp*theta_tmp' + lambda*eye(size(X,2));
%   B = Theta'*Y(i,:)';
%   [L,U,P] = lu(A);
%   x = L*U\(P*B);
%   %x = A\B;
%   X(i,:) = x';
% end
% 
% for i=1:n
%    idx=find(R(:,i)==1);%column 1 entries that == 1
%    x_tmp=X(idx,:);
%    A = zeros(size(Theta,2));%right
%    A = A + x_tmp'*x_tmp + lambda*eye(size(X,2));%
%    B = X'*Y(:,i);
%     [L,U,P] = lu(A);
%      x = L*U\(P*B);
%     Theta(i,:) = x;
%     %Theta(i,:) = A\B;
% end

for i=1:m
  idx=find(R(i,:)==1);
  theta_tmp=(Theta(idx,:));
  A = theta_tmp'*theta_tmp + lambda*eye(size(X,2));
  B = Theta'*Y(i,:)';
  %[L,U,P] = lu(A);
  %x = L*U\(P*B);
  x = A\B;
  X(i,:) = x';
end

for i=1:n
   idx=find(R(:,i)==1);%column 1 entries that == 1
   x_tmp=X(idx,:);
   A = zeros(size(Theta,2));%right
   A = A + x_tmp'*x_tmp + lambda*eye(size(X,2));%
   B = X'*Y(:,i);
    %[L,U,P] = lu(A);
     %x = L*U\(P*B);
    %Theta(i,:) = x;
    Theta(i,:) = A\B;
end


end
