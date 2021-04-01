function [x,w]=gauss_hermite(n)
% Computes Gauss-Hermite quadrature nodes and weights
% - n : number of nodes
% - x : n-by-1 vector of Gauss-Hermite quadrature nodes
% - w : n-by-1 vector of Gauss-Hermite quadrature weights
%
% (c) Aliaksandr Zaretski, 2020

x=sort(roots(flip(a(n))));
w=2^(n-1)*factorial(n)*sqrt(pi)/n^2./H(n-1,x).^2;

end


% Computes n-th order derivative of f(x)=exp(-x^2)
function y=f_pr(n,x)
if n<0
    y=0;
elseif n==0
    y=exp(-x.^2);
else
    y=-2*(n-1)*f_pr(n-2,x)-2*x.*f_pr(n-1,x);
end
end


% Evaluates physicist's Hermite polynomial
function y=H(n,x)
y=(-1)^n*exp(x.^2).*f_pr(n,x);
end


% Computes coefficients of physicist's Hermite polynomial
% H_n(x)=a_{n,0}+a_{n,1}x^1+...+a_{n,n}x^n
function y=a(n)
if n<0
    error('n cannot be negative')
elseif n==0
    y=1;
elseif n==1
    y=[0 2];
else
    y=NaN(1,n+1);
    b=[a(n-1) 0 0];
    y(1)=-b(2);
    for i=2:(n+1)
        y(i)=2*b(i-1)-i*b(i+1);
    end
end
end