%gets the sum of the photon flux and the absorption coefficient for front
%illumination, NOTE: you need to be in the main folder and add this
%function to that path for it to work (dont be in the PEC folder)
function [sumFluxTimesA, absorbedFluxCurrentLimit] = getSumFluxTimesAFront(activelayerthickness)
[PECabsorption, ~, ~, PECabsCoeff, PECactiveLayer, lambdaPEC, maxJsc]  = TransferMatrixPEC(activelayerthickness);
%note: absCoeff is in cm-1

spectrumData = load('spectrum.mat');  %mA.cm-2.nm-1
flux = (spectrumData.fluxLambda).';
fluxRange = (spectrumData.lambda).';
fluxInterp = interp1(fluxRange, flux, lambdaPEC, 'linear', 'extrap');
q = 1.602e-19; %C
photonFlux = fluxInterp./q; %in cm-2.s-1.nm-1

% absorbedFluxmAcm2 = sum(fluxInterp.*absorption(activeLayer, :), 2); %in mA/cm2
% absorbedFlux = absorbedFluxmAcm2*10; %A/m2
absorbedFluxCurrentLimit = maxJsc*10; %A/m2

photonFluxAbsorbed = photonFlux.*PECabsorption(PECactiveLayer, :);
fluxTimesA = (photonFluxAbsorbed.*10000).*(PECabsCoeff(PECactiveLayer,:).*100); 
sumFluxTimesA = sum(fluxTimesA, 2);
end





