%PECJV is in terms of V vs RHE and A/m2
function PECJV = generatePECJV10STH(PECThickness)
%PEC Electrical Model (Assume front illuminated)
deviceName = 'PEC_hematite_front_illuminated'

%--------------------USER INPUT PARAMETERS------------------------------
%Options
showHankinsPlots = 1; %shows the figure that verifies that the model can replicate the same results given in the Hankins et al. paper
plotIndJV = 1; %whether or not you want to plot the JV curves of each individual component of the PEC
showInterpPlots = 1; %shows the interpolated values of each JV curve based on J vector
PECActiveLayerThickness = PECThickness; %in nm

%--------------------Overall Model Parameters---------------------------
vStepSize = 0.01; %step size for voltage that each JV curve for each component will be evaluated on
V = 0:vStepSize:2; %voltage range that each JV curve for each component will be evaluated on
pH = 13.9; %6.7;

deltaU = 0.6; %0.42; %this is the additional thermodynamic potential provided by the external bias
%to provide for the PEC system in question

%----------------------Cathode Parameters-------------------------------
jc0 = 3.1; %AHmmm %exchange current density for the hydrogen evolution reaction (A.m-2)
alpha_c = 0.5; %symmetry factor/ charge transfer coefficient for the reaction
R = 8.314; %universal gas constant (J.mol-1.K-1)
T = 298.15; %temperature (K)
F = 96485; %Faraday's constant (C.mol-1)
Eceq = 0; %pH-dependent equilibrium potential for HER (Assume RHE)
jc_cond = 1.8e6; %conductivity (A.V-1.m-1) %AHmmm
t_c = 0.1/1000; %thickness of cathode in m (need to get this info)

%------------------Photoanode Parameters-----------------------------------
%Phi_surf = ones(size(V)); %surface recombination factor = 1
A = 2.61e-3; %fitting coefficient A %AHmmm
B = 13.3; %fitting coefficient B %AHmmm
Ufb = 0.16; %AHmmm 0.18; %-0.297; %Flat band potential for the anode (V_SHE or V_RHE, make sure to specify in the above code)
e = 1.60217e-19; %elementary charge (C)
n0 = 2.73e27; %3.01e25; %1e24;  %electron donor density (m-3) %AHmmm-Mott-Schottky, Wes
eps0 = 8.85e-12; %Permittivity of free space (A.s.V-1.m-1)
epsR = 38.2; %AHmmm 20.1; %68;  %relative permittivity of the photoanode
ja0 = 5.56e-3; %AHmmm %3.8e-3;  %anode exchange current density A.m-2 

%%%%%% ba = 0.295; %Tafel coefficient from AHmmm
%%%%%% ba =  (alpha_a*96485)/(R*T); %Tafel coefficient
alpha_a = 0.001; % 0.214; % 7.5789e-3; %0.18;0.214;  % ; <--from AHmmm 
Eaeq = 1.23; %pH-dependent equilibrium potential for OER in RHE

%-------------USER INPUT PARAMETERS END---------------------------------


%------------------------CALCULATIONS-----------------------------------
%--------------------Cathode Calculations-------------------------------

%HER Reaction overpotential (use Tafel equation relating echem rxn rate with overpotential)
V_her = V;
%V_her = [-2:vStepSize:min(V), V]; % cathode potential - and extension of V
%to -2 %Brian -- Why is this used and not just V?
eta_c = V_her - Eceq; %overpotential for HER reaction, this is a difference so this is in absolute voltage. 
bc = (-alpha_c*F)/(R*T);
jHER = jc0.*exp(bc.*eta_c);

%Ohmic overpotential
jc_ohm = jc_cond.*(V_her/(t_c)); %this is in absolute voltage as well

%--------------------PhotoAnode Calculations-----------------------------
%Photocurrent - Gaertner-Butler equation
eta_ph = V-Ufb; %"overpotential of the photocurrent"
vMinusUfbPowerHalf = real((V-Ufb).^0.5);
Phi_surf = (A.*exp(B.*(eta_ph)))./(1+(A.*exp(B.*(eta_ph)))); %determine phi_surf utilising fitting method in Hankins et al.    
%Surface recombination is a function of bandbending
%sumFluxTimesA is the attenuation corrected flux density that is actually
%absorbed

[sumFluxTimesA, J_limit] = getSumFluxTimesAFront(PECActiveLayerThickness);
ja_withabs = (((2*e*(sumFluxTimesA^2)*epsR*eps0*(1/n0))^0.5).*vMinusUfbPowerHalf).*Phi_surf; %Core G-B equation with surf recombination
ja_ph = ja_withabs./(1+(ja_withabs./J_limit)); %G-B modified to account for limitation when all absorbed photons are contributing to the photcurrent density


%dark current -- must make sure the electrode potentials are in the same
%Voltage to add them up in parallel and make ja_dark a function of band
%bending - Use Butler-volmer equation in low overpotential approximation
Eaeq_SHE = Eaeq - (0.0592*pH);
eta_a = V - Eaeq_SHE ;
%eta_a = eta_ph + Ufb - Eaeq_SHE ;
ba =  (-alpha_a*F)/(R*T); %Tafel coefficient
ja_dark = ja0.*exp(ba.*(eta_a));

%add the currents together cause they are wired in parallel
ja = ja_dark + ja_ph;


%---------------------Electrolyte Calculations---------------------------
je_cond = 20.6; %conductivity of the electrolyte A.V-1.m-1, AHmmm 1M NaOH, pH13.9
t_e = 2e-3; %thickness of the electrolyte layer in m

je_ohm = je_cond.*(V./t_e);


% %-------------------Substrate Calculations--------------------------------
jm_cond = 0.1*100; %conductivity of the substrate A.V-1.m-1 %arbitrary? Ti AHmmm 1.8e6
t_m = 100e-9; %thickness of the substrate layer in m

jm_ohm = jm_cond.*(V./t_m);


%-----------------Connecting the Components in Series--------------------
%find the component with the limiting current to use that for vector
%maxJ = min([max(je_ohm), max(jm_ohm), max(ja), max(jc_ohm), max(jHER)]); %get the limiting current from all JV curves to determine J vector with
%minJ = max([min(abs(je_ohm)), min(abs(jm_ohm)), min(abs(ja)), min(abs(jc_ohm)), min(abs(jHER))]);
maxJ = 100;
minJ = 3.187;
jMapVect = minJ:((maxJ-minJ)/1000):maxJ; %create J vector

[jc_ohm_unique, jc_ind] = unique(jc_ohm); %voltage contribution for cathode ohmic overpotential
jc_ohm_V = interp1(jc_ohm_unique, V(jc_ind), jMapVect, 'linear', 'extrap');

[jm_ohm_unique, jm_ind] = unique(jm_ohm); %voltage contribution for substrate overpotential
jm_ohm_V = interp1(jm_ohm_unique, V(jm_ind), jMapVect, 'linear', 'extrap');

[je_ohm_unique, je_ind] = unique(je_ohm); %voltage contribution for electrolyte overpotential
je_ohm_V = interp1(je_ohm_unique, V(je_ind), jMapVect, 'linear', 'extrap');

[jHER_unique, jHER_ind] = unique(jHER); %voltage contribution for cathodic reaction overpotential
jHER_V = interp1(jHER_unique, eta_c(jHER_ind), jMapVect, 'linear', 'extrap');

[ja_unique, ja_ind] = unique(ja); %voltage contribution for anodic reaction overpotential
ja_V = interp1(ja_unique, eta_ph(ja_ind), jMapVect, 'linear', 'extrap');

%add the curves together
PECV = deltaU + ja_V + jc_ohm_V + (1.*je_ohm_V) + abs(jHER_V) + jm_ohm_V;
PECJV = [PECV; jMapVect];
%-------------------END OF CALCULATIONS---------------------------------


%------------------PLOTTING OF DATA-------------------------------------
if showHankinsPlots == true
    figure
%     plot(ja_uncorrected(1,:), ja_uncorrected(2,:), 'DisplayName', 'Gartner-Butler No Correction', 'LineWidth',2);
%     hold on
%     plot(V, ja_recomb, 'DisplayName', 'GB corrected for Recombination', 'LineWidth',2);
%     hold on
    plot(V, ja_withabs, 'DisplayName', 'GB corrected for Recomb and abs', 'LineWidth',2);
    hold on
    plot(V,ja_ph, 'DisplayName', 'GB corrected for recomb, abs, and with photon flux limitation', 'LineWidth',2);
    grid on
    legend
    xlabel('Applied Electrode Potential (V_S_H_E)')
    ylabel('j A/m2')
    title('J_p_h with different factors accounted (as shown in Hankins et al.)')
    
    annadarkfile = fopen('darkcurrentanna.txt');
    annaDCdata = textscan(annadarkfile, '%f %f', 'HeaderLines', 0);
    annaDC = [annaDCdata{1}'; annaDCdata{2}'];
    fclose(annadarkfile);
    
    annaPCfile = fopen('photocurrentanna.txt');
    annaPCdata = textscan(annaPCfile, '%f %f', 'HeaderLines', 0);
    annaPC = [annaPCdata{1}'; annaPCdata{2}'];
    fclose(annaPCfile);
   
    figure
    plot(V, ja_dark, 'DisplayName', append('Modelled Anodic Dark Current with \alpha_a = ', string(alpha_a)),'LineWidth', 1.5)
    hold on
    plot(V, ja_ph, 'DisplayName', 'Modelled Anodic Photocurrent','LineWidth', 1.5)
    hold on
    plot(annaDC(1,:), annaDC(2,:), 'DisplayName', 'Experimentally Recorded Dark Current (Hankins et al.)', 'LineWidth', 3, 'LineStyle', ':')
    hold on
    plot(annaPC(1,:), annaPC(2,:), 'DisplayName', 'Experimentally Recorded  Photocurrent (Hankins et al.)', 'LineWidth', 3, 'LineStyle', ':')
    legend
    grid on
    axis([0 max(V) 0 max(ja_ph)])
    xlabel('Electrode Potential (E_a) (V_S_H_E)')
    ylabel('Current Density, j (A/m2)')
    axis([0 1.4 0 40])
    title('Comparison Between Modelled J_a_,_d_a_r_k and J_p_h and experimental J_a_,_d_a_r_k and J_p_h from Hankins et al. ')
    
end

if plotIndJV == true
    figure
    plot(eta_c, jHER,'LineWidth',2);
    grid on
    ylabel('Current Density, | j |, (A/m2)')
    xlabel('Overpotential, \eta_c (V)')
    title('Cathode HER Overpotential')
    
    figure
    plot(V, jc_ohm, 'LineWidth',2);
    grid on 
    ylabel('Current Density, J, (A/m2)')
    xlabel('Overpotential, j.d_c/\sigma_c (V)')
    title('Cathode Ohmic Overpotential')
    
    figure
    plot(V, je_ohm,'LineWidth',2);
    grid on 
    ylabel('Current Density, J, (A/m2)')
    xlabel('Overpotential, j.d_e/\sigma_e (V)')
    title('Electrolyte Ohmic Overpotential')
%     
%     figure
%     plot(V, jm_ohm,'LineWidth',2);
%     grid on 
%     ylabel('Current Density, J, (A/m2)')
%     xlabel('Overpotential, j.d_m/\sigma_m (V)')
%     title('Membrane Ohmic Overpotential')
    
    figure
    plot(eta_ph, ja_dark, 'DisplayName', 'Anodic Dark Current','LineWidth', 1.5)
    hold on
    plot(eta_ph, ja_ph, 'DisplayName', 'Anodic Photocurrent','LineWidth', 1.5)
    hold on
    plot(eta_ph, ja, 'DisplayName', 'Anodic Current Density (J_p_h + j_d_a_r_k)', 'LineWidth', 1.5)
    legend
    grid on
    axis([0 max(V) 0 max(ja_ph)])
    xlabel('\eta_p_h (E_a - U_f_b) (V)')
    ylabel('Current Density, j (A/m2)')
    
%     BiVO4file = fopen('BiVO4_DC.txt');
%     BiVO4DCData = textscan(BiVO4file, '%f %f', 'HeaderLines', 0);
%     BiVO4DC = [BiVO4DCData{1}.'; BiVO4DCData{2}.'];
%     fclose(BiVO4file);
    
    figure
    plot(V, ja_dark, 'DisplayName', 'Anodic Dark Current','LineWidth', 1.5)
    hold on
    plot(V, ja_ph, 'DisplayName', 'Anodic Photocurrent','LineWidth', 1.5)
    hold on
    plot(V, ja, 'DisplayName', 'Anodic Current Density (J_p_h + j_d_a_r_k)', 'LineWidth', 1.5)
%     hold on
%     plot(BiVO4DC(1,:), BiVO4DC(2,:), 'DisplayName', 'Durant et al. Dark Current Data', 'LineWidth', 1.5)
    legend
    grid on
    axis([0 max(V) 0 max(ja_ph)])
    xlabel('Electrode Potential (E_a) (V_S_H_E)')
    ylabel('Current Density, j (A/m2)')
    axis([0 1.5 0 50])

      
end


if showInterpPlots == true
    %plots to compare the interpolation
    figure
%     plot(V, jc_ohm, 'DisplayName', 'Jc non interpolated', 'LineWidth',2);
%     hold on
    plot(jc_ohm_V, jMapVect, 'DisplayName', 'Jc interpolated', 'LineWidth',2);
    ylabel('Current Density, J, (A/m2)')
    xlabel('Overpotential, j.d_c/\sigma_c  (V)')
    grid on
    legend
    title('Cathode Ohmic Overpotential Interpolated')
    
    figure
%     plot(V, je_ohm, 'DisplayName', 'Je non interpolated', 'LineWidth',2);
%     hold on
    plot(je_ohm_V, jMapVect, 'DisplayName', 'Je interpolated', 'LineWidth',2);
    ylabel('Current Density, J, (A/m2)')
    xlabel('Overpotential, j.d_e/\sigma_e (V)')
    grid on
    legend
    title('Electrolyte Overpotential Interpolated')
    
    figure
%     plot(V, jm_ohm, 'DisplayName', 'Jm non interpolated', 'LineWidth',2);
%     hold on
    plot(jm_ohm_V, jMapVect, 'DisplayName', 'Jm interpolated', 'LineWidth',2);
    ylabel('Current Density, J, (A/m2)')
    xlabel('Overpotential, j.d_m/\sigma_m (V)')
    grid on
    legend
    title('Membrane Overpotential Interpolated')
    
    figure
%     plot(eta_c, jHER, 'DisplayName', 'jHER non interpolated', 'LineWidth',2);
%     hold on
    plot(jHER_V, jMapVect, 'DisplayName', 'jHER interpolated', 'LineWidth',2);
    ylabel('Current Density, J, (A/m2)')
    xlabel('HER Overpotential, \eta_c (V)')
    grid on
    legend
    title('HER JV Curve Interpolated')
    
    figure
%     plot(eta_ph, ja, 'DisplayName', 'Anode Current Density non-interped', 'LineWidth',2);
%     hold on
    plot(ja_V, jMapVect, 'DisplayName', 'Anode Current Density interpolated', 'LineWidth',2);
    ylabel('Current Density, J, (A/m2)')
    xlabel('Overpotential, \Delta\phi_S_C (V)')
    grid on
    legend
    title('Anode JV Curve Interpolated')
    
end

%determine the point where current density for the system is 8.14mA/cm2
STHVoltage = interp1(PECJV(2,:), PECJV(1,:), 81.301, 'linear', 'extrap');
PhotoanodeVoltage = interp1(jMapVect, ja_V, 81.301, 'linear', 'extrap');

figure
plot(PECJV(1,:), PECJV(2,:), 'DisplayName', 'PEC System JV', 'LineWidth',1.5);
hold on
plot(jc_ohm_V, jMapVect, 'DisplayName', 'Cathode Ohmic Overpotential JV', 'LineWidth',1.5);
hold on
plot(jHER_V, jMapVect, 'DisplayName', 'HER Overpotential JV', 'LineWidth',1.5);
hold on
plot(je_ohm_V, jMapVect, 'DisplayName', 'Electrolyte Overpotential JV', 'LineWidth',1.5);
hold on
plot(jm_ohm_V, jMapVect, 'DisplayName', 'Substrate Ohmic Overpotential JV', 'LineWidth',1.5);
hold on
plot(ja_V, jMapVect, 'DisplayName', 'Photoanode Current Density (J_p_h + J_d_a_r_k)', 'LineWidth', 1.5);

ylin1 = yline(81.301, '--k', { '81.301 A/m^2'},'LineWidth', 1.5);
ylin1.Annotation.LegendInformation.IconDisplayStyle = 'off';
ylin1.LabelHorizontalAlignment = 'left';
xlin1 = xline(STHVoltage, '--k', {'STH 10% Voltage =', append(string(STHVoltage), ' V')},'LineWidth', 1.5);
xlin1.Annotation.LegendInformation.IconDisplayStyle = 'off';
xlin1.LabelOrientation = 'horizontal';
xlin1.LabelVerticalAlignment = 'bottom';

xlin2 = xline(PhotoanodeVoltage, '--k', {'Photoanode Voltage =', append(string(PhotoanodeVoltage), ' V')},'LineWidth', 1.5);
xlin2.Annotation.LegendInformation.IconDisplayStyle = 'off';
xlin2.LabelOrientation = 'horizontal';
xlin2.LabelVerticalAlignment = 'bottom';

ylabel('Current Density, | j | (A/m2)')
xlabel('Overpotential (V) / Photoelectrochemical Cell Potential (V)')
title('J-V Curve for the PEC System and its Components')
%axis([-0.2 2 minJ maxJ]); 
axis([-0.2 2 minJ 100]); 
legend 
grid on
%-------------------END OF PLOTTING DATA------------------------------

% save(append(fileName, '_data'));
% simulationFigures = findobj('Type','Figure');
% savefig(simulationFigures, append( fileName, '_figures.fig'));

end 
