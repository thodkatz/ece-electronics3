% AEM = 9282
AEM = 82;

% Specs
fprintf("Specs\n")
CL = (2 + AEM*0.01);   % pf
SR = 18 + 0.01*AEM;    % V/uS
Vdd = 1.8 + 0.003*AEM; % V
Vss = -Vdd;            % V
GB = 7 + 0.01*AEM;     % MHz
A = 20 + 0.01*AEM;     % dB
P = 50 + 0.01*AEM;     % mW
fprintf("CL: %f\nSR: %f\nVdd: %f\nVss: %f\n", CL, SR, Vdd, Vss);
fprintf("GB: %f\nA: %f\nP: %f\n", GB, A, P);
fprintf("\n")

% Vin
Vin_max  = 100e-3;
Vin_min = -100e-3;

% constants
Eox = 3.45e-11;

% NMOS
Ln   = 0.05; % parameter035.pdf
Vtn  = 0.768;
%Vtn_max = ?
%Vtn_min = ?
Uoxn = 591.7 * 1e-4; % m^2/V*s
Toxn = 2.12e-8;
Coxn = Eox/Toxn;
Kn   = Uoxn * Coxn; % manual
Kn_spice = 9.693e-5; % spice parameter
fprintf("NMOS\n");
fprintf("Cox: %f\nKn: %f (manual)\nKn: %f (spice param)\n",Coxn,Kn,Kn_spice);

fprintf("\n")

% PMOS
Lp   = 0.15; % parameter035.pdf 
Vtp  = -0.9056;
%Vtp_max = ?
%Vtp_min = ?
Uoxp = 180.2 * 1e-4;
Toxp = 2.12e-8;
Coxp = Eox/Toxp;
Kp   = Uoxp * Coxp; % manual
Kp_spice = 2.9352e-5; % spice parameter
fprintf("PMOS\n")
fprintf("Cox: %f\nKp: %f (manual)\nKp: %f (spice param)\n", Coxp, Kp, Kp_spice);
fprintf("\n")

% Step 1
L = 1e-6; %um
fprintf("L: %fm, %fum\n\n", L, L*1e6);

% Step 2
Cc = 0.22*CL; %pF
Cc = ceil(Cc);

% Step 3
I5 = SR*Cc * 1e-6; %Ampere
fprintf("I5: %fA, %fuA\n", I5, I5 * 1e6)

% Step 4
S3 = I5 / ((Kp) * (Vdd - Vin_max - abs(Vtp) + Vtn)^2);
fprintf("S3: %f\n", S3)
if S3 < 1
    S3 = 1;
end
S4 = S3;
fprintf("S3 and S4: %f\n\n", S3)

% Step 5
I3 = I5/2;
gm3  = sqrt(2 * Kp * S3 * I3);
W3 = S3 * L;
Cgs3 = 2 * 0.667 * W3 * L * Coxp;
p3   = gm3/(2*Cgs3); 
fprintf("P3: %fMHz\n", (p3/(2*pi))*1e-6)
assert(p3/(2*pi) > 10*GB);

% Step 6
gm1 = GB*1e6 * 2*pi * Cc * 1e-12; %uS
gm2 = gm1; % CMOS
S1  = gm1^2/(Kn * I5);
S1  = ceil(S1);
S2  = S1;

% Step 7
b1   = S1*Kn;
Vds5 = Vin_min - Vss - sqrt(I5/b1) - Vtn;
if Vds5 < 0.1
    fprintf("Vds: %f\nPick a bigger value for Vds5\n", Vds5);
end
S5 = (2*I5)/(Kp*(Vds5)^2);
S5 = ceil(S5); % try to unceil this, then S6 != S7

% Step 8
gm6 = 2.2*gm2*(CL/Cc);
if (gm6 < 10*gm1)
   gm6 = 10*gm1;
end
I4 = I3;
gm4 = sqrt(2*Kp*S4*I4);
S6  = S4*(gm6/gm4);
S6 = ceil(S6);

% Do we have vout max spec? No
% Step 9
%Vds6 = Vdd - Vout
%temp = gm6/(Kp*Vds6)
%S6   = max(S6, temp)

I6 = (gm6^2)/(2*Kp*S6);

% Step 10
S7 = (I6/I5)*S5;
S7 = ceil(S7);
b7 = S7 * Kn;
I7 = I6;
%Vds7 = sqrt((2*I7)/b7) - Vtp;
%fprintf("Vds7: %f and Vout: %f\n", Vds7, Vout);
%assert(Vds7 < Vout);

% Step 11
Av = (2*gm2*gm6)/(I5*(Ln + Lp) * I6*(Ln + Lp));
Pd = (I5 + I6)*(Vdd + abs(Vss));

fprintf("Gain: %f\nGoal: %f\n", Av, A);
if (mag2db(Av) < A)
    fprintf("You need to increase the gain\n");
elseif (Pd *1e3 > P)
    fprintf("Power dissipation: %f. Goal: %f\nYou need to decrease the power dissipation\n", Pd, P);
end

%um
W1 = S1 * L * 1e6;
W2 = S2 * L * 1e6;
W3 = S3 * L * 1e6;
W4 = W3;
W5 = S5 * L * 1e6;
W6 = S6 * L * 1e6;
W7 = S7 * L * 1e6;
W8 = W5;

W = [W1 W2 W3 W4 W5 W6 W7 W8]
Iref = I5;
I7   = I6;
I = [I3*1e6 I4*1e6 I5*1e6 I6*1e6 I7*1e6 Iref*1e6] %uA


