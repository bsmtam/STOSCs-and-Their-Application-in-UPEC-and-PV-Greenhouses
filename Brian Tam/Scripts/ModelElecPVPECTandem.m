
% Models the JV curves for the PV and PEC cell and then connects them in
% series electrically to get an operating point
close all
clear
clc
fileName = 'nil';
PVactivelayerthickness = 0; %nm
PECactivelayerthikness = 0; %nm

isFrontIlluminated = 1;

%Generate PV JV curve, uncomment/comment the relevant PVJV scripts you want
%to use
if isFrontIlluminated == true
    PVJV = generatePVJVFront(PECactivelayerthikness, PVactivelayerthickness);
    %PVJV = generatePVJVFrontSTH10(PECactivelayerthikness);
elseif isFrontIlluminated == false
    PVJV = generatePVJVBack;
end


%Generate PEC Curve
PECJV = generatePECJV3(isFrontIlluminated, PECactivelayerthikness, PVactivelayerthickness);
%PECJV = generatePECJV10STH(isFrontIlluminated, PECactivelayerthikness);

%change current density values to mA/cm2 from A/m2
PECJV(2,:) = PECJV(2,:).*0.1;


%find operating point (connect them in tandem)
%define range for new J vector to ensure no crazy extrapolations to
%minimise error
maxJ = min([max(PVJV(2,:)); max(PECJV(2,:)) ]); %largest J value in the new J vector is defined by JV curve with the smaller upper limit
minJ = max([min(PVJV(2,:)); min(PECJV(2,:))]); %smallest J value in the new J vector is defined by the JV curve with the largest lower limit
jMapVector = linspace(minJ, maxJ, 1000);


%interpolate them on the same J vector
[PV_J_unique, PVJ_ind] = unique(PVJV(2,:));
PV_V = interp1(PV_J_unique, PVJV(1,PVJ_ind), jMapVector, 'linear', 'extrap');
%PV_V = interp1(PVJV(2,:), PVJV(1,:), jMapVector, 'linear', 'extrap');
PEC_V = interp1(PECJV(2,:), PECJV(1,:), jMapVector, 'linear', 'extrap');

PVPECTandem_JV = [(PV_V-PEC_V); jMapVector];

jOpt = interp1(PVPECTandem_JV(1,:),PVPECTandem_JV(2,:), 0, 'linear');
vOpt = interp1(PECJV(2,:), PECJV(1,:), jOpt, 'linear');

figure
plot(PVJV(1,:), PVJV(2,:), 'DisplayName', 'PV Cell', 'LineWidth', 1.5);
hold on
plot(PECJV(1,:), PECJV(2,:), 'DisplayName', 'PEC Cell', 'LineWidth', 1.5);
legend
ylin1 = yline(jOpt, '--k', {'Operating Current Density =', append(string(jOpt), ' mA/cm^2')},'LineWidth', 1.5);
ylin1.Annotation.LegendInformation.IconDisplayStyle = 'off';
ylin1.LabelHorizontalAlignment = 'left';
xlin1 = xline(vOpt, '--k', {'Operating Voltage =', append(string(vOpt), ' V')},'LineWidth', 1.5);
xlin1.Annotation.LegendInformation.IconDisplayStyle = 'off';
xlin1.LabelOrientation = 'horizontal';
xlin1.LabelVerticalAlignment = 'bottom';
grid on
axis([0 max(PVJV(1,:)) 0 max(PVJV(2,:))])
title('J-V Curves for the PEC system and PV Cell')
xlabel('Potential (V)')
ylabel('J (mA/cm2)')


figure
plot(PVPECTandem_JV(1,:), PVPECTandem_JV(2,:), 'LineWidth', 1.5);
xlin = xline(0, '--k', 'V = 0', 'LineWidth', 1.5);
xlin.LabelOrientation = 'horizontal';
yline(jOpt, '-r', {'Operating Current Density =', append(string(jOpt), ' mA/cm^2')}, 'LineWidth', 1.5)
title('PV-PEC tandem JV Curve');
xlabel('V')
ylabel('J (mA/cm2)')
grid on



%Calculate STH
etaF = 1; %faradaic efficiency
Pin = 0.1; %W/cm2
STH = ((1.23*jOpt/1000*etaF)/Pin)*100;

fprintf('STH is %.2f%%', STH);

if isFrontIlluminated == 1
    save(append(fileName, '_PVPECTandemFront_data'));
    simulationFigures = findobj('Type','Figure');
    savefig(simulationFigures, append( fileName, ' _PVPECTandemFront_figures.fig'));
    
elseif isFrontIlluminated == 0
    save(append(fileName, '_PVPECTandemBack_data'));
    simulationFigures = findobj('Type','Figure');
    savefig(simulationFigures, append( fileName, ' _PVPECTandemBack_figures.fig'));
end
