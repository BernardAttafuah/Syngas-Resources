%%Assumptions
%1.Non-Isothermal Packed Bed Reactor Model
%2.Negligible pressure drop across the reactor
%3.No shaft work
%4.Adiabatic Reactor Operation

domain = [0 0.5];
initial_temperature = 2900;
initial_conversion = 0;
IC = [initial_temperature ,initial_conversion];

[w,X] = ode45(@nopressuredrop,domain,IC);
subplot(2,1,1)
plot(w,X(:,1),'b')
hold on
xlabel('Weight')
ylabel('Temperature')
title('Temperature Profile')
grid on
hold off
hold on
subplot(2,1,2)
plot(w,X(:,2),'k')
hold on
xlabel('Weight')
ylabel('Conversion')
title('Conversion Profile')
grid on
%Testing to see if a point will lie on the graph without considering
%non-isothermal
WeightTest = 0.211;
ConversionTest = 0.8988;
hold on
%plot(WeightTest,ConversionTest,'g^')
hold off
grid on
hold off
grid on

function [nonisothermal] = nopressuredrop(w,X)
%Constants Parameters
T_Reference = 298;
FAo=1000;
FBo=1500;
To=1900;
Po=30*10^5;
CpCO=32.455;
CpH2O=39.93;
CpH2=29.915;
CpCH4=67.184;
standard_enthalpy_reaction= 205900;
R=8.314;
yAo=FAo/(FBo+FAo);
CAo= (yAo*Po)/(R*To);
delta_Cp= (3*CpH2)+CpCO-CpH2O-CpCH4;
Phi_B=FBo/FAo;
epsilon = yAo*((3+1)-(1+1));
CA=(CAo*(1-X(2))*To)/(X(1)*(1+(epsilon*X(2))));
CB=(CAo*(Phi_B-X(2))*To)/(X(1)*(1+(epsilon*X(2))));
k=(3*10^5)*exp(-15107/X(1));
Rate=k*CA*CB;
H_rxn_T=standard_enthalpy_reaction+ delta_Cp*(X(1)-T_Reference);
nonisothermal =[(-H_rxn_T*Rate)/(FAo*(delta_Cp*X(2)+(CpCH4+Phi_B*CpH2O)));(Rate/FAo)];
end




