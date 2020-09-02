close all
clear
clc

load('GHProcessingData.mat');

[a, b] = size(GHResultsProcessing);

%b = 10;
linewidth = 2;
%create Jsc vs active layer thickness plot for emmott
figure
for i = 1:1:b
    
%     if i ~= 12 && i~= 14  && i~= 15 && i~= 17 && i~= 21 && i~= 23
%         
%         plot(thicknessRange, [GHResultsProcessing(i).Jsc], 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth)
%         hold on
%         
%     end

switch i
    case 1
        plot(thicknessRange, [GHResultsProcessing(i).Jsc], 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', 'r')
        hold on
    case 2
        plot(thicknessRange, [GHResultsProcessing(i).Jsc], 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', 'r', 'LineStyle', '--')
        hold on
    case 3
        plot(thicknessRange, [GHResultsProcessing(i).Jsc], 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', 'b')
        hold on
    case 4
        plot(thicknessRange, [GHResultsProcessing(i).Jsc], 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', 'b', 'LineStyle', '--')
        hold on
    case 5
        plot(thicknessRange, [GHResultsProcessing(i).Jsc], 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', [0.4940 0.1840 0.5560])
        hold on
    case 6
        plot(thicknessRange, [GHResultsProcessing(i).Jsc], 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', [0.4940 0.1840 0.5560], 'LineStyle', '--')
        hold on
    case 7
        plot(thicknessRange, [GHResultsProcessing(i).Jsc], 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', [0.9290 0.6940 0.1250])
        hold on
    case 8
        plot(thicknessRange, [GHResultsProcessing(i).Jsc], 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', [0.9290 0.6940 0.1250], 'LineStyle', '--')
        hold on
    case 9
        plot(thicknessRange, [GHResultsProcessing(i).Jsc], 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', 'k')
        hold on
    case 10
        plot(thicknessRange, [GHResultsProcessing(i).Jsc], 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', 'k', 'LineStyle', '--')
        hold on
end 
    
    
end
legend
grid on
xlabel('Active Layer Thickness (nm)')
ylabel('Jsc (mA/cm2)')
title('Jsc Against Active Layer Thickness for Devices in the Emmott et al. study')

%Create G Vs layer thickness for emmott
figure
for i = 1:1:b
    
%     if i ~= 12 && i~= 14  && i~= 15 && i~= 17 && i~= 21 && i~= 23
%         
%         plot(thicknessRange, [GHResultsProcessing(i).G].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth)
%         hold on
%         
%     end

switch i
    case 1
        plot(thicknessRange, [GHResultsProcessing(i).G].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', 'r')
        hold on
    case 2
        plot(thicknessRange, [GHResultsProcessing(i).G].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', 'r', 'LineStyle', '--')
        hold on
    case 3
        plot(thicknessRange, [GHResultsProcessing(i).G].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', 'b')
        hold on
    case 4
        plot(thicknessRange, [GHResultsProcessing(i).G].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', 'b', 'LineStyle', '--')
        hold on
    case 5
        plot(thicknessRange, [GHResultsProcessing(i).G].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', [0.4940 0.1840 0.5560])
        hold on
    case 6
        plot(thicknessRange, [GHResultsProcessing(i).G].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', [0.4940 0.1840 0.5560], 'LineStyle', '--')
        hold on
    case 7
        plot(thicknessRange, [GHResultsProcessing(i).G].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', [0.9290 0.6940 0.1250])
        hold on
    case 8
        plot(thicknessRange, [GHResultsProcessing(i).G].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', [0.9290 0.6940 0.1250], 'LineStyle', '--')
        hold on
    case 9
        plot(thicknessRange, [GHResultsProcessing(i).G].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', 'k')
        hold on
    case 10
        plot(thicknessRange, [GHResultsProcessing(i).G].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', 'k', 'LineStyle', '--')
        hold on
end 
    
    
end
legend
grid on
xlabel('Active Layer Thickness (nm)')
ylabel('Growth Factor (%)')
title('Growth Factor Against Active Layer Thickness for Devices in the Emmott et al. study')

%create G vs PCE emmott

figure
for i = 1:1:b
    
%     if i ~= 12 && i~= 14  && i~= 15 && i~= 17 && i~= 21 && i~= 23
%         
%         plot(thicknessRange, [GHResultsProcessing(i).PCE].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth)
%         hold on
%         
%     end

switch i
    case 1
        plot([GHResultsProcessing(i).PCE].*100,[GHResultsProcessing(i).G].*100,  'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', 'r')
        hold on
    case 2
        plot([GHResultsProcessing(i).PCE].*100,[GHResultsProcessing(i).G].*100,  'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', 'r', 'LineStyle', '--')
        hold on
    case 3
        plot([GHResultsProcessing(i).PCE].*100,[GHResultsProcessing(i).G].*100,  'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', 'b')
        hold on
    case 4
        plot([GHResultsProcessing(i).PCE].*100,[GHResultsProcessing(i).G].*100,  'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', 'b', 'LineStyle', '--')
        hold on
    case 5
        plot([GHResultsProcessing(i).PCE].*100,[GHResultsProcessing(i).G].*100,  'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', [0.4940 0.1840 0.5560])
        hold on
    case 6
        plot( [GHResultsProcessing(i).PCE].*100,[GHResultsProcessing(i).G].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', [0.4940 0.1840 0.5560], 'LineStyle', '--')
        hold on
    case 7
        plot( [GHResultsProcessing(i).PCE].*100,[GHResultsProcessing(i).G].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', [0.9290 0.6940 0.1250])
        hold on
    case 8
        plot([GHResultsProcessing(i).PCE].*100,[GHResultsProcessing(i).G].*100,  'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', [0.9290 0.6940 0.1250], 'LineStyle', '--')
        hold on
    case 9
        plot( [GHResultsProcessing(i).PCE].*100,[GHResultsProcessing(i).G].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', 'k')
        hold on
    case 10
        plot( [GHResultsProcessing(i).PCE].*100,[GHResultsProcessing(i).G].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', 'k', 'LineStyle', '--')
        hold on
end 
    
    
end
legend
grid on
xlabel('PCE (%)')
ylabel('Growth Factor (%)')
title('Growth Factor Against PCE for Devices in the Emmott et al. study')
xlin3 = xline(2.1, '--k', {'Best Case PCE Required = 2.1%'},'LineWidth', 1.5);
xlin3.Annotation.LegendInformation.IconDisplayStyle = 'off';
xlin3.LabelOrientation = 'aligned';
xlin3.LabelVerticalAlignment = 'bottom';
xlin4 = xline(1.31, '--k', {'Very Low Cost OPV PCE Required = 1.31%'},'LineWidth', 1.5);
xlin4.Annotation.LegendInformation.IconDisplayStyle = 'off';
xlin4.LabelOrientation = 'aligned';
xlin4.LabelVerticalAlignment = 'bottom';

% Plot new device data

figure
for i = 1:1:b

switch i
    case 11
        plot(thicknessRange, [GHResultsProcessing(i).Jsc], 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', 'r')
        hold on
    case 13
        plot(thicknessRange, [GHResultsProcessing(i).Jsc], 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', 'r', 'LineStyle', '--')
        hold on
    case 16
        plot(thicknessRange, [GHResultsProcessing(i).Jsc], 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', 'b')
        hold on
    case 18
        plot(thicknessRange, [GHResultsProcessing(i).Jsc], 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', 'b', 'LineStyle', '--')
        hold on
    case 19
        plot(thicknessRange, [GHResultsProcessing(i).Jsc], 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', [0.4940 0.1840 0.5560])
        hold on
    case 20
        plot(thicknessRange, [GHResultsProcessing(i).Jsc], 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', [0.9290 0.6940 0.1250])
        hold on
    case 22
        plot(thicknessRange, [GHResultsProcessing(i).Jsc], 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', 'k')
        hold on
    case 24
        plot(thicknessRange, [GHResultsProcessing(i).Jsc], 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', [0 0.5 0])
        hold on
    case 5
        plot(thicknessRange, [GHResultsProcessing(i).Jsc], 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', 3, 'Color', [0 0.75 0.75], 'LineStyle', '-.')
        hold on
end 
    
    
end
legend
grid on
xlabel('Active Layer Thickness (nm)')
ylabel('Jsc (mA/cm2)')
title('Jsc Against Active Layer Thickness for More Recent Devices')


figure
for i = 1:1:b

switch i
    case 11
        plot(thicknessRange, [GHResultsProcessing(i).G].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', 'r')
        hold on
    case 13
        plot(thicknessRange, [GHResultsProcessing(i).G].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', 'r', 'LineStyle', '--')
        hold on
    case 16
        plot(thicknessRange, [GHResultsProcessing(i).G].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', 'b')
        hold on
    case 18
        plot(thicknessRange, [GHResultsProcessing(i).G].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', 'b', 'LineStyle', '--')
        hold on
    case 19
        plot(thicknessRange, [GHResultsProcessing(i).G].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', [0.4940 0.1840 0.5560])
        hold on
    case 20
        plot(thicknessRange, [GHResultsProcessing(i).G].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', [0.9290 0.6940 0.1250])
        hold on
    case 22
        plot(thicknessRange, [GHResultsProcessing(i).G].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', 'k')
        hold on
    case 24
        plot(thicknessRange, [GHResultsProcessing(i).G].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', [0 0.5 0])
        hold on
    case 5
        plot(thicknessRange, [GHResultsProcessing(i).G].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', 3, 'Color', [0 0.75 0.75], 'LineStyle', '-.')
        hold on
    
end 
    
    
end
legend
grid on
xlabel('Active Layer Thickness (nm)')
ylabel('Growth Factor (%)')
title('Growth Factor Against Active Layer Thickness for More Recent Devices')

figure
for i = 1:1:b

switch i
    case 11
        plot([GHResultsProcessing(i).PCE].*100, [GHResultsProcessing(i).G].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', 'r')
        hold on
    case 13
        plot([GHResultsProcessing(i).PCE].*100, [GHResultsProcessing(i).G].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', 'r', 'LineStyle', '--')
        hold on
    case 16
        plot([GHResultsProcessing(i).PCE].*100, [GHResultsProcessing(i).G].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', 'b')
        hold on
    case 18
        plot([GHResultsProcessing(i).PCE].*100, [GHResultsProcessing(i).G].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', 'b', 'LineStyle', '--')
        hold on
    case 19
        plot([GHResultsProcessing(i).PCE].*100, [GHResultsProcessing(i).G].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', [0.4940 0.1840 0.5560])
        hold on
    case 20
        plot([GHResultsProcessing(i).PCE].*100, [GHResultsProcessing(i).G].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', [0.9290 0.6940 0.1250])
        hold on
    case 22
        plot([GHResultsProcessing(i).PCE].*100, [GHResultsProcessing(i).G].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', 'k')
        hold on
    case 24
        plot([GHResultsProcessing(i).PCE].*100, [GHResultsProcessing(i).G].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', linewidth, 'Color', [0 0.5 0])
        hold on
    case 5
        plot([GHResultsProcessing(i).PCE].*100, [GHResultsProcessing(i).G].*100, 'DisplayName', GHResultsProcessing(i).deviceName, 'LineWidth', 3, 'Color', [0 0.75 0.75], 'LineStyle', '-.')
        hold on
    
end 
    
    
end
legend
grid on
xlabel('PCE (%)')
ylabel('Growth Factor (%)')
title('Growth Factor Against PCE for More Recent Devices')
xlin1 = xline(2.1, '--k', {'Best Case PCE Required = 2.1%'},'LineWidth', 1.5);
xlin1.Annotation.LegendInformation.IconDisplayStyle = 'off';
xlin1.LabelOrientation = 'aligned';
xlin1.LabelVerticalAlignment = 'bottom';
xlin2 = xline(1.31, '--k', {'Very Low Cost OPV PCE Required = 1.31%'},'LineWidth', 1.5);
xlin2.Annotation.LegendInformation.IconDisplayStyle = 'off';
xlin2.LabelOrientation = 'aligned';
xlin2.LabelVerticalAlignment = 'bottom';

%plot solar insolation graph
cSiEff = 0.2;
CIGSeff = 0.12;
perfRatio = 0.8;

csIInsoNE = [];
CIGsInsoNE = [];

csIInsoSE = [];
CIGsInsoSE = [];

NEInsolation = 1300*10000/1000; %MWh/ha
SEInsolation = 2200*10000/1000; %MWh/ha

for coverage = 0:0.05:1
    csIInsoNE = [csIInsoNE, NEInsolation*perfRatio*coverage*cSiEff];
    CIGsInsoNE = [CIGsInsoNE, NEInsolation*perfRatio*coverage*CIGSeff];
end

for coverage = 0:0.05:1
    csIInsoSE  = [csIInsoSE , SEInsolation*perfRatio*coverage*cSiEff];
    CIGsInsoSE = [CIGsInsoSE, SEInsolation*perfRatio*coverage*CIGSeff];
end

OpaqueGFactor = (1 - [0:0.05:1])*0.91;

CuSCNY6PCE = [GHResultsProcessing(13).PCE];
CuSCNY6NEInso = CuSCNY6PCE.*NEInsolation.*perfRatio;
CuSCNY6SEInso = CuSCNY6PCE.*SEInsolation.*perfRatio;
CuSCNY6G = [GHResultsProcessing(13).G.*0.91];
% CuSCNY6G = CuSCNY6G + 0.07;

PMDPP3TPCE = [GHResultsProcessing(5).PCE];
PMDPP3TPCENEInso = PMDPP3TPCE.*NEInsolation.*perfRatio;
PMDPP3TPCESEInso = PMDPP3TPCE.*SEInsolation.*perfRatio;
PMDPP3TG = [GHResultsProcessing(5).G.*0.91];

PETGNE = ones(1, length(csIInsoNE)).*91;
PETGSE = ones(1, length(csIInsoSE)).*91;

figure
plot(csIInsoNE, OpaqueGFactor.*100, 'DisplayName', 'opaque c-Si', 'LineWidth', 2, 'LineStyle', '--')
hold on
plot(CIGsInsoNE,OpaqueGFactor.*100, 'DisplayName', 'opaque CIGS', 'LineWidth', 2, 'LineStyle', '--' )
hold on
plot(CuSCNY6NEInso, CuSCNY6G.*100, 'DisplayName', 'CuSCN:Y6', 'LineWidth', 2)
hold on
plot(PMDPP3TPCENEInso, PMDPP3TG.*100, 'DisplayName', 'PMDPP3T:PC60BM', 'LineWidth', 2)
hold on
plot(csIInsoNE, PETGNE, 'DisplayName', 'Clear glass or PET', 'LineWidth', 3, 'LineStyle', ':')
xlabel('Electricity Yield (MWh/ha/year)')
ylabel('Growth Factor (%)')
title('Growth Factor Against Annual Electricity Yield for Northern Europe for different PV greenhouse tech ')
legend
grid on

figure
plot(csIInsoSE, OpaqueGFactor.*100, 'DisplayName', 'opaque c-Si', 'LineWidth', 2, 'LineStyle', '--')
hold on
plot(CIGsInsoSE,OpaqueGFactor.*100, 'DisplayName', 'opaque CIGS', 'LineWidth', 2, 'LineStyle', '--' )
hold on
plot(CuSCNY6SEInso, CuSCNY6G.*100, 'DisplayName', 'CuSCN:Y6', 'LineWidth', 2)
hold on
plot(PMDPP3TPCESEInso, PMDPP3TG.*100, 'DisplayName', 'PMDPP3T:PC60BM', 'LineWidth', 2)
hold on
plot(csIInsoSE, PETGSE, 'DisplayName', 'Clear glass or PET', 'LineWidth', 3, 'LineStyle', ':')
xlabel('Electricity Yield (MWh/ha/year)')
ylabel('Growth Factor (%)')
title('Growth Factor Against Annual Electricity Yield for Southern Europe for different PV greenhouse tech ')
legend
grid on

DataFigures = findobj('Type','Figure');
savefig(DataFigures, 'GreenhouseStudy_figures.fig');