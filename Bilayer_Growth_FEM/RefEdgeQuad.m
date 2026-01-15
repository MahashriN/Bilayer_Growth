function [w, p] = RefEdgeQuad(order)
% Returns quadrature weights w and points p for 1D interval [0,1]
    switch order
        case 1  % Midpoint rule
            w = 1.0;
            p = 0.5;
        case 2  % 2-point Gauss quadrature
            p = [0.21132486540519, 0.78867513459481];
            w = [0.5, 0.5];
        case 3  % 3-point Gauss quadrature
            p = [0.11270166537926, 0.5, 0.88729833462074];
            w = [5/18, 8/18, 5/18];
        otherwise
            error('Edge quadrature order not implemented');
    end
end
