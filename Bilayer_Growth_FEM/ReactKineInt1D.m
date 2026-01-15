function F = ReactKineInt1D(p1, edges1, p2, edges2, u1, v1, u2, v2, eta, f, ord)
% ReactKineInt1D computes ∫ f(u(x)) * φ_i(x) dx over 1D domain
% Inputs:
%   p     - 1 x N node coordinates
%   edges - Ne x 2 connectivity of nodes (each row = [i j])
%   u     - N x 1 nodal values
%   f     - function handle for reaction kinetics (e.g., f = @(u) a - u)
%   ord   - quadrature order

N = length(p1);          % Number of nodes
Ne1 = size(edges1, 1);    % Number of elements (edges)
F = zeros(N, 1);        % Global load vector

% Quadrature rule and basis functions
[iw, ip] = RefEdgeQuad(ord);             % weights and points on [0,1]
[basis_ip,~] = basis_linear_1D(ip);          % 2 x n_quad_pts

% Loop over each element (edge)
for k = 1:Ne1
    nodes1 = edges1(k, :);
    nodes2 = edges2(k, :);
    x1 = p1(nodes1);                    % Coordinates of the two nodes
    x2 = p2(nodes2);
    h = abs(x1(2) - x1(1));               % Element length

    u1_local = u1(nodes1);             % Local nodal values
    v1_local = v1(nodes1);
    u2_local = u2(nodes2);
    v2_local = v2(nodes2);
    u1_ip = u1_local' * basis_ip;         % u at quadrature points (1 x n_quad)
    v1_ip = v1_local' * basis_ip;
    u2_ip = u2_local' * basis_ip;
    v2_ip = v2_local' * basis_ip;
    
    f_vals = f(u1_ip,v1_ip,u2_ip,v2_ip,eta);
    f_vals(1);
if any(isnan(u1_ip)) || any(isnan(u2_ip)) || any(isnan(v1_ip)) || any(isnan(v2_ip))
    error('NaN detected in u_ip values at element %d', k);
end

if any(isnan(f_vals))
    disp('Inputs to f:');
    disp([u1_ip; v1_ip; u2_ip; v2_ip]);
    error('f returned NaN at element %d', k);
end
    % Compute local load vector
    local_load = zeros(2, 1);
    for i = 1:2
        local_load(i) = sum(f_vals .* basis_ip(i,:) .* iw) * h;
    end

    % Assemble into global load vector
    F(nodes1) = F(nodes1) + local_load;
end
end
