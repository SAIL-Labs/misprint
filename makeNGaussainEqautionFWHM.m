function eqnStr=makeNGaussainEqautionFWHM(N)
eqnStr='a1*exp(- (x-b1).^2 / (2*(fwhm1/(2*sqrt(2*log(2))))^2))';  % stddev

for i=2:N
   eqnStr=[eqnStr ' + a' num2str(i) '*exp(- (x-b' num2str(i) ').^2 / (2*(fwhm' num2str(i) '/(2*sqrt(2*log(2))))^2))'];
end
eqnStr=[eqnStr '+ d'];



