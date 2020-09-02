fileName = 'PC61BM';

CSVTrans = readmatrix(append(fileName, '_Transmission.csv'), 'HeaderLines', 0);
CSVTransModel = readmatrix(append(fileName, '_Model_Transmission.csv'), 'HeaderLines', 0);

figure
plot(CSVTrans(:,1), CSVTrans(:,2), 'DisplayName', 'Experimental Transmission','LineWidth', 4)
hold on
plot(CSVTransModel(:,1), CSVTransModel(:,2), 'DisplayName', 'Modelled Transmission Without Ag','LineWidth', 2)
hold on
grid on
xlabel('Wavelength (nm)')
ylabel('Transmission (%)')
title('Transmission Profiles of Experimental Transmission and Modelled Transmission Without Ag')
legend


CSVTransInterp = [CSVTransModel(:,1),interp1(CSVTrans(:,1), CSVTrans(:,2), CSVTransModel(:,1),'linear','extrap')]
plot(CSVTransInterp(:,1), CSVTransInterp(:,2), 'DisplayName', 'Interpolated Exp. Values','LineWidth', 2)

AgFactor = (CSVTransInterp(:,2)./CSVTransModel(:,2)).';

hold on
plot(CSVTransModel(:,1), CSVTransModel(:,2).*AgFactor.','DisplayName', 'Model with Ag Factor','LineWidth', 1, 'LineStyle',':' )


save('PC61BM_d1_AgFactor.mat', 'AgFactor')
simulationFigures = findobj('Type','Figure');
savefig(simulationFigures, append( fileName, '_figures.fig'));