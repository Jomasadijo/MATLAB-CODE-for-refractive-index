clear;

%% Input values
n1 = 1.74742;
n2 = 1.82866;
n3 = 2.1085;
n4 = 1.9612;
n5 = 1.46;
n6 = 1;

n9 = 1.394;
n10 = 0.8395 + 6.1276i;

% Layer thicknesses
d1 = 60e-9;
d2 = 40e-9;
d3 = 10e-9;
d4 = 150e-9;
d5 = 1.1e-3;

d9 = 1e-9;
d10 = 100e-9;

%% Constants
c_const = 3e8;  % Speed of light

% Wavelength (can adjust for different cases)
lam0 = 520e-9;

% Refractive indices for the materials (adjustable)
eps1 = (n1)^2;
eps2 = (n2)^2;
eps3 = (n3)^2;
eps4 = (n4)^2;
eps5 = (n5)^2;
eps6 = (n6)^2;

eps9 = (n9)^2;
eps10 = (n10)^2;

%% Wavenumber and k values
u = 0.00001:0.0001:2;
k1 = 2*pi/lam0 * sqrt(eps1);
k1x = u .* k1;
k1z = sqrt(k1^2 - k1x.^2);

%% Reflection coefficients (for TM and TE)
l1 = -i*sqrt(eps1/eps1 - u.^2);  
l2 = -i*sqrt(eps2/eps1 - u.^2);  
l3 = -i*sqrt(eps3/eps1 - u.^2);  
l4 = -i*sqrt(eps4/eps1 - u.^2);  
l5 = -i*sqrt(eps5/eps1 - u.^2);  
l6 = -i*sqrt(eps6/eps1 - u.^2);
deltaph2 = k1 * 2 * l2 * d2;  
deltaph3 = k1 * 2 * l3 * d3;  
deltaph4 = k1 * 2 * l4 * d4; 

% TM and TE coefficients
r1p = (l2*eps1 - l1*eps2) ./ (l2*eps1 + l1*eps2);  
r1s = (l2 - l1) ./ (l2 + l1);  
t1s = 2 * l2 ./ (l1 + l2);
t1p = 2 * sqrt(eps1) .* sqrt(eps2) .* l2 ./ (l2 * eps1 + l1 * eps2);

r2p = (l3*eps2 - l2*eps3) ./ (l3*eps2 + l2*eps3);   
r2s = (l3 - l2) ./ (l3 + l2);  
t2s = 2 * l3 ./ (l2 + l3);
t2p = 2 * sqrt(eps2) .* sqrt(eps3) .* l3 ./ (l2 * eps3 + l3 * eps2);

% Continue calculating higher-order reflection coefficients
r3p = (l4*eps3 - l3*eps4) ./ (l4*eps3 + l3*eps4);  
r3s = (l4 - l3) ./ (l4 + l3);  
t3s = 2 * l4 ./ (l4 + l3);
t3p = 2 * sqrt(eps3) .* sqrt(eps4) .* l4 ./ (l3 * eps4 + l4 * eps3);

r4p = (l5*eps4 - l4*eps5) ./ (l5*eps4 + l4*eps5);  
r4s = (l5 - l4) ./ (l5 + l4);   
t4s = 2 * l5 ./ (l5 + l4);
t4p = 2 * sqrt(eps4) .* sqrt(eps5) .* l5 ./ (l4 * eps5 + l5 * eps4);

r5p = (l6*eps5 - l5*eps6) ./ (l6*eps5 + l5*eps6);  
Rp_air = abs(r5p).^2;
r5s = (l6 - l5) ./ (l6 + l5);   
Rs_air = abs(r5s).^2;
t5s = 2 * l6 ./ (l5 + l6);  
Ts_air = abs(t5s).^2;
t5p = 2 * sqrt(eps5) .* sqrt(eps6) .* l6 ./ (l5 * eps6 + l6 * eps5);  
Tp_air = abs(t5p).^2;

% Define rp_f for the layer interface (for example, layer between eps6 and air or some other material)
rp_f = (l6*eps6 - l5*eps5) ./ (l6*eps6 + l5*eps5); % Example definition

%% Reflection and transmission at interfaces
rpp1 = (r3p + r4p .* exp(-deltaph4)) ./ (1 + r3p .* r4p .* exp(-deltaph4));    
rss1 = (r3s + r4s .* exp(-deltaph4)) ./ (1 + r3s .* r4s .* exp(-deltaph4));

rpp2 = (r2p + rpp1 .* exp(-deltaph3)) ./ (1 + r2p .* rpp1 .* exp(-deltaph3));  
rss2 = (r2s + rss1 .* exp(-deltaph3)) ./ (1 + r2s .* rss1 .* exp(-deltaph3));

rp_z = (r1p + rpp2 .* exp(-deltaph2)) ./ (1 + r1p .* rpp2 .* exp(-deltaph2));  
rs_z = (r1s + rss2 .* exp(-deltaph2)) ./ (1 + r1s .* rss2 .* exp(-deltaph2));

Rp_zong = abs(rp_z).^2;
Rs_zong = abs(rs_z).^2;

%% Transmission coefficients
tpp1 = (t3p .* t4p .* exp(-deltaph4 ./ 2)) ./ (1 + r3p .* r4p .* exp(-deltaph4));     
tss1 = (t3s .* t4s .* exp(-deltaph4 ./ 2)) ./ (1 + r3s .* r4s .* exp(-deltaph4));

tpp2 = (t2p .* tpp1 .* exp(-deltaph3 ./ 2)) ./ (1 + r2p .* rpp1 .* exp(-deltaph3)); 
tss2 = (t2s .* tss1 .* exp(-deltaph3 ./ 2)) ./ (1 + r2s .* rss1 .* exp(-deltaph3));

tTM_z = (t1p .* tpp2 .* exp(-deltaph2 ./ 2)) ./ (1 + r1p .* rpp2 .* exp(-deltaph2)); 
tTE_z = (t1s .* tss2 .* exp(-deltaph2 ./ 2)) ./ (1 + r1s .* rss2 .* exp(-deltaph2));

% More coefficients for final transmission terms
beta_z = k1 * 2 * l1 * 0;   %(出光方向)
beta_f = k1 * 2 * l1 * d1;  %(阴极方向)

% Correcting the final expression for KTMv
KTMv = 3 / 4 * real(((1 - rp_z .* exp(-beta_z)) .* (1 - rp_f .* exp(-beta_f))) ...
    ./ (1 - rp_z .* rp_f .* exp(-(beta_z + beta_f))));  % Correct parentheses

% Define the URL of the image
imageURL = 'https://qph.cf2.quoracdn.net/main-qimg-1d05bfbbffd4dc62cd2c5d3a3795a1d8';

% Read the image data from the URL
imageData = webread(imageURL);
imshow(imageData);
title('Fraction of Power vs. Ag-EML Distance');

%% Output values (Formatted)

disp('Reflection and Transmission Coefficients:');
disp('----------------------------------------');

disp('Reflection Coefficients (Rp_zong):');
disp(['Rp_zong: ', num2str(Rp_zong(1))]);

disp('Reflection Coefficients (Rs_zong):');
disp(['Rs_zong: ', num2str(Rs_zong(1))]);

disp('Transmission Coefficients (tTM_z):');
disp(['tTM_z: ', num2str(real(tTM_z(1)))]);
disp(['tTM_z (Imaginary part): ', num2str(imag(tTM_z(1)))]);

disp('Transmission Coefficients (tTE_z):');
disp(['tTE_z: ', num2str(real(tTE_z(1)))]);
disp(['tTE_z (Imaginary part): ', num2str(imag(tTE_z(1)))]);

disp('Final KTMv value:');
disp(['KTMv: ', num2str(real(KTMv(1)))]);

%% Scatter Plot with Curves (Blue and Red)
figure;
hold on;

% Scatter plot
scatter(u, Rp_zong, 'b', 'filled'); % Blue scatter for Rp_zong
scatter(u, Rs_zong, 'r', 'filled'); % Red scatter for Rs_zong

% Curve for Rp_zong
plot(u, Rp_zong, 'b-', 'LineWidth', 2); % Blue line for Rp_zong

% Curve for Rs_zong
plot(u, Rs_zong, 'r-', 'LineWidth', 2); % Red line for Rs_zong

title('Reflection Coefficients (Rp_zong & Rs_zong)');
xlabel('u');
ylabel('Coefficient Value');
legend({'Rp_zong', 'Rs_zong'}, 'Location', 'best');
grid on;
hold off;
