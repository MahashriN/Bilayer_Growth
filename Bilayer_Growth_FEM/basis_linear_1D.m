function [basis, grad_basis] = basis_linear_1D(x)
% Evaluates linear basis functions and their gradients on [0,1]
x = x(:)';  % Ensure row vector for consistent output
n = length(x);

basis = [1 - x; x];        % φ₁ = 1 - x, φ₂ = x

grad_basis = repmat([-1; 1], 1, n);  % dφ₁/dξ = -1, dφ₂/dξ = 1
end
