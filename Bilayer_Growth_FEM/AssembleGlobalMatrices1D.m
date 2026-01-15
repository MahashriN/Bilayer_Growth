function [St, Ma] = AssembleGlobalMatrices1D(p, t, ord)
% Assembles 1D stiffness and mass matrices using linear basis functions

N = length(p);         % Number of nodes
M = size(t, 1);        % Number of elements
St = sparse(N, N);
Ma = sparse(N, N);

% Get quadrature points and weights on reference interval [0,1]
[iw, ip] = RefEdgeQuad(ord);     % 1D quadrature
[basis_ip,grad_basis] = basis_linear_1D(ip);  % shape: 2 x num_quad_pts

% Loop over each 1D element
for k = 1:M
    node_ids = t(k,:);           % Node indices of current element
    x = p(node_ids);              % Node coordinates (2x1)

    h = abs(x(2) - x(1));         % Length of the interval

    % Derivatives of basis functions in real coordinates
    grad_basis_real = grad_basis / h;

    % Local matrices
    Local_M = zeros(2, 2);
    Local_S = zeros(2, 2);

    for q = 1:length(ip)
        phi = basis_ip(:, q);     % Basis at quadrature point
        wq = iw(q);               % Weight at this point

        % Mass matrix contribution
        Local_M = Local_M + (phi * phi') * wq * h;

        % Stiffness matrix contribution (constant for linear basis)
        % So we can compute this outside the loop (optional)
    end

    % Constant stiffness matrix for linear basis
    Local_S = (grad_basis_real * grad_basis_real') * h;

    % Assemble into global matrices
    Ma(node_ids, node_ids) = Ma(node_ids, node_ids) + Local_M;
    St(node_ids, node_ids) = St(node_ids, node_ids) + Local_S;
end
end
