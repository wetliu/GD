function [J, X_grad, Theta_grad] = cofiCostFunc4(X, Theta, Y, R, lambda)
%COFICOSTFUNC Collaborative filtering cost function
%   [J, grad] = COFICOSTFUNC(params, Y, R, num_users, num_movies, ...
%   num_features, lambda) returns the cost and gradient for the
%   collaborative filtering problem.

            
% You need to return the following values correctly
J = 0;
X_grad = zeros(size(X));
Theta_grad = zeros(size(Theta));
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
J = 0;
J=1/2*sum(sum(R.*((X*Theta'-Y).^2)))+lambda/2*sum(sum(Theta.^2))+lambda/2*sum(sum(X.^2));
k = size(X,2);
for i=1:size(X,1)
  idx=find(R(i,:)==1);
  for p = 1:size(idx,2)
      resi = Y(i,idx(p)) - X(i,:)*Theta(idx(p),:)';
      %X(i,:) = (resi) * Theta(idx(p),:)./(lambda+Theta(idx(p),:).^2);
    %for j = 1:k
    X_grad(i,:) = 2*( resi)*Theta(idx(p),:); %/ (lambda+Theta(idx(p),j)^2);%Theta(idx(p),j)*X(i,j)
      Theta_grad(idx(p),:) = 2*(resi)*X(i,:); %/ (lambda+X(i,j)^2);%Theta(idx(p),j)*X(i,j)
      %X_grad(i,j) = ( resi)*Theta(idx(p),j); %/ (lambda+Theta(idx(p),j)^2);%Theta(idx(p),j)*X(i,j)
      %Theta_grad(idx(p),j) = (resi)*X(i,j); %/ (lambda+X(i,j)^2);%Theta(idx(p),j)*X(i,j)
  end
end


end
