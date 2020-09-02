%gets the sum of the photon flux and the absorption coefficient for front
%illumination, NOTE: you need to be in the main folder and add this
%function to that path for it to work (dont be in the PEC folder)
function [sumFluxTimesA, absorbedFluxCurrentLimit] = getSumFluxTimesABack(PECactivelayerthickness, PVactivelayerthickness)
fileName = 'nil';

[PECabsorption, ~, ~, PECabsCoeff, PECactiveLayer, lambdaPEC, maxJsc]  = TransferMatrixPEC(PECactivelayerthickness);
[~, ~,PVtransmission,~, ~, lambdaPV]= TransferMatrixPV(PVactivelayerthickness); 

PVTransInterp = interp1(lambdaPV, PVtransmission, lambdaPEC, 'linear', 'extrap');
PVTransInterp(PVTransInterp < 0) = 0;

totalabsorption = PVTransInterp.*PECabsorption(PECactiveLayer,:);

figure
plot(lambdaPEC, PECabsorption(PECactiveLayer,:).*100, 'DisplayName', 'PEC absorption under AM1.5', 'LineWidth', 2)
hold on
plot(lambdaPV, PVtransmission.*100, 'DisplayName', 'PV transmission under AM1.5', 'LineWidth', 2)
hold on
plot(lambdaPV, totalabsorption.*100, 'DisplayName', 'Modified PEC Absorption', 'LineWidth', 2)
xlabel('Wavelength (nm)')
ylabel('%')
grid on
legend
title('Absorption and Transmission Spectra for the PEC and the PV, with Modified PEC Absorption spectrum')

spectrumData = load('spectrum.mat');  %mA.cm-2.nm-1
flux = (spectrumData.fluxLambda).';
fluxRange = (spectrumData.lambda).';
fluxInterp = interp1(fluxRange, flux, lambdaPEC, 'linear', 'extrap');
q = 1.602e-19; %C
photonFlux = fluxInterp./q; %in cm-2.s-1.nm-1

%absorbedFluxmAcm2 = sum(fluxInterp.*totalabsorption, 2); %in mA/cm2
absorbedFluxmAcm2 = trapz(lambdaPEC,fluxInterp.*totalabsorption); %in mA/cm2
absorbedFluxCurrentLimit = absorbedFluxmAcm2*10; %A/m2


photonFluxAbsorbed = photonFlux.*totalabsorption;
fluxTimesA = (photonFluxAbsorbed.*10000).*(PECabsCoeff(PECactiveLayer,:).*100);
sumFluxTimesA = sum(fluxTimesA, 2);

save(append(fileName, '_getSumFluxTimesABack_data'));
    simulationFigures = findobj('Type','Figure');
    savefig(simulationFigures, append( fileName, ' _getSumFluxTimesABack_figures.fig'));
end





