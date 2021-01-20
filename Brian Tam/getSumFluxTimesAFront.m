%gets the sum of the photon flux multiplied with absorption coefficient for front
%illumination
function [sumFluxTimesA, absorbedFluxCurrentLimit] = getSumFluxTimesAFront(activelayerthickness)
[PECabsorption, ~, ~, PECabsCoeff, PECactiveLayer, lambdaPEC, maxJsc]  = TransferMatrixPEC(activelayerthickness);
%note: absCoeff is in cm-1

spectrumData = load('spectrum.mat');  %mA.cm-2.nm-1 solar spectrum
flux = (spectrumData.fluxLambda).';
fluxRange = (spectrumData.lambda).'; %nm
%fluxInterp = interp1(fluxRange, flux, lambdaPEC, 'linear', 'extrap'); %mA.cm-2.nm-1 
% q = 1.602e-19; %CT
% photonFlux = fluxInterp./q; %in 0.001 cm-2.s-1.nm-1

%%%%% -Brian - Assume flux starts as mW.cm-2.nm-1 = mJ.cm-2.s-1.nm-1
%%%%% instead of %mA.cm-2.nm-1 
%Divide by plank constant (Js), divide by c (m/s), multiply by wavelength
%of light (m)
photonFlux_preinterp = flux./6.626e-34/3e8.*fluxRange/1e9; %in 0.001 cm-2.s-1.nm-1
photonFlux = interp1(fluxRange, photonFlux_preinterp, lambdaPEC, 'linear', 'extrap'); %mA.cm-2.nm-1 

% absorbedFluxmAcm2 = sum(fluxInterp.*absorption(activeLayer, :), 2); %in mA/cm2
% absorbedFlux = absorbedFluxmAcm2*10; %A/m2
absorbedFluxCurrentLimit = maxJsc*10; %A/m2 (was in mA/cm^2) originally

% Absorption coefficient in cm^-1 (JAP Vol 86 p.487 Eq 23)
photonFluxAbsorbed = photonFlux.*PECabsorption(PECactiveLayer, :); %in 0.001 cm-2.s-1.nm-1
fluxTimesA = (photonFluxAbsorbed.*10000).*(PECabsCoeff(PECactiveLayer,:).*100); %in m-2.s-1.nm-1.m-1
%%%% Brian - change in units here there is *10000 to convert from cm-2 to
%%%% m-2 but also a *0.001 to convert from mA to A
%fluxTimesA = (photonFluxAbsorbed.*10000).*(PECabsCoeff(PECactiveLayer,:).*100);
sumFluxTimesA = sum(fluxTimesA, 2); %in m-2.s-1.nm-1.m-1
end





