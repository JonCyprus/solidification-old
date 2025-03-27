function [ dfx  ] = nablax(f,delta)

% three point stencil gradient
%dfx =  (circshift( f, [0, -1] ) - circshift( f, [ 0, 1]))/(2 * delta);

% five-point stencil
dfx = (-circshift(f, [0, -2]) + 8 * circshift(f, [0, -1]) ...
           - 8 * circshift(f, [0, 1]) + circshift(f, [0, 2])) / (12 * delta);
end

