weightspan = linspace(0,0.5,10);
Xo = 0;
[v,X] = ode45(@nopressure_isothermal,weightspan,Xo);
%graph
plot(v,X,'b')
xlabel('Weight in kg')
ylabel('Conversion')
title('A graph of Conversion against weight  -(Isothermal with Negligible Pressure drop)')
grid on
legend('Conversion Profile')

%Testing to see if a point will lie on the graph
hold on

WeightTest = 0.211;
ConversionTest = 0.8988;
hold on
plot(WeightTest,ConversionTest,'r^')
hold off
grid on


function dXdv = nopressure_isothermal(v,X)
FAo=1000;
FBo=1500;
To=1000;
Po=30*10^5;
R=8.314;
yAo=FAo/(FBo+FAo)
epsilon = yAo*((3+1)-(1+1))
CAo = (yAo*Po)/(R*To)
Phi_B=FBo/FAo
CA=(CAo*(1-X))/(1+(epsilon*X));
CB=(CAo*(Phi_B-X))/(1+(epsilon*X));
k=1.49;
Rate=k*CA*CB;
dXdv = (Rate/FAo);
end