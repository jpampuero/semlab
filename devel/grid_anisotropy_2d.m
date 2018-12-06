th=linspace(-pi,pi,1000);
for p=[2 4 6],
  polar(th, cos(th).^(2*p+2) +sin(th).^(2*p+2) )
  hold all
end
legend('2','4','6')
hold off
