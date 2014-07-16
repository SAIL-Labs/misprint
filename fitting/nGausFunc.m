function out=nGausFunc(x,xData,N)
% x is format  [amps means sigmas(:) offset], as col vect.

if nargin==2
    N=(length(x0)-1)/3;
end
x(end)=0;

out=bsxfun(@times,x(1:N), exp(- bsxfun(@rdivide, [bsxfun(@minus,repmat(xData,[1,N]),x(N+1:N*2)).^2],[2*x(N*2+1:N*3).^2]) )) + x(end);