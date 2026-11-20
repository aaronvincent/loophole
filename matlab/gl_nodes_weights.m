function [x, w] = gl_nodes_weights(n)
% gl_nodes_weights  n-point Gauss-Legendre nodes/weights on [-1,1].
% Golub-Welsch: eigen-decomposition of the symmetric tridiagonal Jacobi
% matrix for Legendre polynomials.  int_{-1}^{1} f dx ~ sum(w .* f(x)).
persistent cache
if isempty(cache), cache = containers.Map('KeyType','double','ValueType','any'); end
if isKey(cache, n)
    xw = cache(n); x = xw{1}; w = xw{2}; return;
end
k    = (1:n-1);
beta = k ./ sqrt(4*k.^2 - 1);
J    = diag(beta,1) + diag(beta,-1);
[V,D] = eig(J);
[x,idx] = sort(diag(D));
w = 2 * (V(1,idx).^2).';
x = x(:); w = w(:);
cache(n) = {x, w};
end
