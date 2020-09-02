
function PVJV = generatePVJVBack
%-Model utilises parameters to plot out JV curve for a back illuminated
%solar cell in a PV-PEC tandem (light hits PV first) this essentially would
%just use the Jsc and Voc reported in literature to create a curve 

%------------------Parameters from Literature at STC-----------------
%get these 2 values from the device literature
fileName = 'nil';
fileName2 = 'nil';
Jsc = 1; %Jsc = Jph (assumed) (mA/cm2)
Voc = 1; %Voc (V)

%perhaps get these values from digitising curves? idk.
%get shunt and series resistance from a JV curve
JVfile = fopen(append(fileName, '_JV.txt'));
JVdata = textscan(JVfile, '%f %f', 'HeaderLines', 0);
expJV = [JVdata{1}'; JVdata{2}']; % experimental JV curve with 1st row in V, 2nd row in mA/cm2
fclose(JVfile);

%interpolate the curve to have as many points as possible
R_volt_range = linspace(0, Voc, 1000);
expJV_interp = [R_volt_range; interp1(expJV(1,:), expJV(2,:), R_volt_range, 'linear', 'extrap')];

%get gradient vector of the JV curve
gradExpJV = gradient(expJV_interp(2,:), expJV_interp(1,:));

%shunt resistance is estimated at the slope near the Jsc
Rsh = 1/abs(gradExpJV(1)); %shunt resistance (milliOhms/cm2)



%series resistance is estimated at the slope near the Voc
Rs = 1/abs(gradExpJV(length(gradExpJV))); %series resistance (milliOhms/cm2)

%--------------------------------------------------------------------

k = 1.3806488e-23; %boltzmann's constant (J/K)
T = 298.15; %temperature in K
n = 1.5; %ideality factor
q = 1.60217662e-19; %elementary charge (C)
voltageStepSize = 0.01;

%Calculated Parameters
Vt = (k*T)/q; %thermal voltage (V)
V  = 0:voltageStepSize:Voc+0.2;
J0 = Jsc/(exp(Voc/(n*Vt))-1);

JV = @(J) Jsc - J0.*(exp((V + J.*Rs)./(n*Vt))-1) - (V + J.*Rs)./Rsh - J;

J_guess = Jsc*ones(size(V));
J_sol = fsolve(JV, J_guess);

PVJV = [V; J_sol];

figure
plot(V, J_sol, 'LineWidth', 2, 'DisplayName', 'Model');
hold on
plot(expJV(1,:), expJV(2,:), 'DisplayName', 'Experimental data from literature', 'LineWidth', 2);
xlabel('Voltage (V)')
ylabel('J (mA/cm2)')
axis([0 (Voc+voltageStepSize) 0 (Jsc+10*voltageStepSize)])
title('JV Curve for PV Cell')
grid on
legend

%compare FF
expPower = [expJV_interp(1,:); expJV_interp(1,:).*expJV_interp(2,:)];
[maxPexp, maxPexInd] = max(expPower(2,:));
expFF = maxPexp/(Jsc*Voc)
expVm = expJV_interp(1, maxPexInd)
expJm = expJV_interp(2, maxPexInd)

modPower = [V; J_sol.*V];
[maxP, maxPInd] = max(modPower(2,:))
modFF = maxP/(Jsc*Voc)
modVm = V(maxPInd)
modJm = J_sol(maxPInd)

save(append(fileName2, '_generatePVJVBack_data'));
simulationFigures = findobj('Type','Figure');
savefig(simulationFigures, append( fileName2, '_generatePVJVBack_figures.fig'));
end

