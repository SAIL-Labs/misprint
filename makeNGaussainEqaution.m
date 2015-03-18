function eqnStr=makeNGaussainEqaution(N)
eqnStr='a1*exp(- (x-b1).^2 / (2*c1^2 ))';  % stddev

for i=2:N
   eqnStr=[eqnStr ' + a' num2str(i) '*exp(- (x-b' num2str(i) ').^2 / (2*c' num2str(i) '.^2 ))'];
end
eqnStr=[eqnStr '+ d'];