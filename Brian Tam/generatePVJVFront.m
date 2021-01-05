%-Model utilises parameters to plot out JV curve for a front illuminated
%solar cell in a PV-PEC tandem (light hits PEC first)

%methodology: 
%1. get EQE from lit and digitise that
%2. get absorption profile from simulating using TFM
%3. get IQE from that
%4. get new absorption profile by multiplying it with the transmission of
%the PEC and get the adjusted EQE
%5. get a new Jsc by integrating the the product of photonflux and EQE over
%lambda

function PVJV = generatePVJVFront(PECthickness, PVThickness)

%------------------Parameters from Literature at STC-----------------
%get these 2 values from the device literature (values for STC conditions)
fileName = 'nil';
fileName2 = 'nil';
Voc_lit = 1; %Voc at AM1.5 STC (given in lit)
PVActiveLayerThickness = PVThickness; %nm
PECActiveLayerThickness = PECthickness; %nm
makingTandem = false;

k = 1.3806488e-23; %boltzmann's constant (J/K)
T = 298.15; %temperature in K
n = 1.5; %ideality factor
q = 1.60217662e-19; %elementary charge (C)

%get shunt and series resistance from a JV curve
JVfile = fopen(append(fileName, '_JV.txt'));
JVdata = textscan(JVfile, '%f %f', 'HeaderLines', 0);
expJV = [JVdata{1}'; JVdata{2}']; % experimental JV curve with 1st row in V, 2nd row in mA/cm2
fclose(JVfile);

%interpolate the curve to have as many points as possible
R_volt_range = linspace(0, Voc_lit, 1000);
expJV_interp = [R_volt_range; interp1(expJV(1,:), expJV(2,:), R_volt_range, 'linear', 'extrap')];

%get gradient vector of the JV curve
gradExpJV = gradient(expJV_interp(2,:), expJV_interp(1,:));

%shunt resistance is estimated at the slope near the Jsc
Rsh = 1/abs(gradExpJV(1)); %shunt resistance (milliOhms/cm2)



%series resistance is estimated at the slope near the Voc
Rs = 1/abs(gradExpJV(length(gradExpJV))); %series resistance (milliOhms/cm2)


%get adjusted Jsc
[PVabsorption, ~,~,~, PVactiveLayer, lambdaPV ]= TransferMatrixPV(PVActiveLayerThickness); 
[~, ~, PEC_transmission, ~, ~, lambdaPEC, ~]  = TransferMatrixPEC(PECActiveLayerThickness);

%make sure the EQE files are in terms of nm and for the EQE make sure
%you're clear whether its in % or in fractions because Flurin's ITIC and Y6
%data is in % so you gotta divide it by 100 for the Jsc later
EQEfile = fopen(append(fileName, '_EQE.txt'));
EQEdata = textscan(EQEfile, '%f %f', 'HeaderLines', 0);
EQE = [EQEdata{1}'; EQEdata{2}']; %first row is lambda, second row is the EQE
fclose(EQEfile);

%get IQE
EQEinterp = interp1(EQE(1,:), EQE(2,:), lambdaPV, 'linear', 'extrap');
EQEinterp(EQEinterp<0) = 0; %make sure the interpolation doesnt go negative
IQE = EQEinterp./PVabsorption(PVactiveLayer, :);
IQE(IQE==inf) = 0;
IQE(IQE<0) = 0;
IQE(isnan(IQE)) = 0;

% figure
% plot(EQE(1,:), EQE(2,:), 'DisplayName', 'not interped')
% hold on
% plot(lambdaPV, EQEinterp, 'DisplayName', 'interped')
% legend
% grid on

%get adjusted absorption profile
PECtransInterp = interp1(lambdaPEC, PEC_transmission, lambdaPV, 'linear', 'extrap');
%PECtransInterp = ones(1, length(lambdaPV)).*0.8; %sanity check
totalAbsorption = PECtransInterp.*PVabsorption(PVactiveLayer, :);

%new EQE
EQEmod = totalAbsorption.*IQE;

figure
plot(lambdaPV, EQEmod,'DisplayName' ,'modified EQE', 'LineWidth', 2);
hold on
plot(EQE(1,:), EQE(2,:), 'DisplayName', 'Unmodified digitised EQE', 'LineWidth', 2);
hold on
plot(lambdaPV, EQEinterp, 'DisplayName', 'Interpolated digitised EQE', 'LineStyle', ':')
hold on
plot(lambdaPEC, totalAbsorption.*100, 'DisplayName', 'Modified Absorption of PV Cell', 'LineWidth', 2)
hold on
plot(lambdaPV, PVabsorption(PVactiveLayer,:).*100, 'DisplayName', 'PV active layer absorption', 'LineWidth', 2)
hold on
plot(lambdaPEC, PEC_transmission.*100, 'DisplayName', 'PEC transmission', 'LineWidth', 2)
hold on
plot(lambdaPEC, PECtransInterp.*100, 'DisplayName', 'Interpolated PEC transmission', 'LineStyle', ':')
legend
xlabel('Wavelength (nm)')
ylabel('%')
title('Absorption/Transmission Profiles and Unmodified/Modified EQE Spectra')
grid on

%get adjusted Jsc
data = load('spectrum.mat'); % in ./Data
fluxLambda=data.fluxLambda; %in mA/cm2.nm1
lambda1=data.lambda;

fluxLambdaInterp = interp1(lambda1, fluxLambda, lambdaPV, 'linear', 'extrap');
Jsc = trapz(lambdaPV,EQEmod.*fluxLambdaInterp)./100; %in mA/cm2

if makingTandem == true
    Jsc = Jsc/2;
end

%get adjusted Voc for light intensity
irradianceAM15 = trapz(lambdaPV, fluxLambdaInterp.*(1240./lambdaPV));
incidentIrradiance = trapz(lambdaPV,PECtransInterp.*fluxLambdaInterp.*(1240./lambdaPV));
concFactor = incidentIrradiance/irradianceAM15;
Voc = Voc_lit + ((n*k*T)/q)*log(concFactor); %adjusted Voc for intensity difference


%------------------------------------------------------------------------


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
plot(V, J_sol, 'LineWidth', 2, 'DisplayName', 'Modelled JV Curve After Accounting for PEC Absorption');
hold on
xlabel('Voltage (V)')
ylabel('J (mA/cm2)')
axis([0 (Voc+voltageStepSize) 0 (Jsc+5*voltageStepSize)])
title(append('JV Curve for PV Cell: ', fileName))
grid on

save(append(fileName2, '_generatePVJVFront_data'));
simulationFigures = findobj('Type','Figure');
savefig(simulationFigures, append( fileName2, '_generatePVJVFront_figures.fig'));


end

