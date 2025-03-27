function [ dfy  ] = nablay(f,delta)

% three point stencil gradient
%dfy =  (circshift( f, [-1, 0] ) - circshift( f, [1, 0]))/(2 * delta);

% five point stencil gradient
dfy = (-circshift(f, [-2,0]) + 8 * circshift(f, [-1, 0]) ...
           - 8 * circshift(f, [1, 0]) + circshift(f, [2, 0])) / (12 * delta);
end

