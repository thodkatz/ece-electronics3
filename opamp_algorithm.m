% AEM = 9282
AEM = 82;

% Specs
CL = 2 + AME*0.01;
SR = 18 + 0.01*AEM;
Vdd = 1.8 + 0.003*AEM;
Vss = -Vdd;
GB = 7 + 0.01*AEM;
A = 20 + 0.01*AEM;
P = 50 + 0.01*AEM;

% Vin, Vout (max)
Vin  = 100e-3;
Vout = 100e-3;

% constants
Eox = 3.45e-11;

% NMOS
Ln   = 0.04;
Vtn  = 0.768;
Uoxn = 591.7;
Toxn = 2.12e-8;
Coxn = Eox/Toxn;
Kn   = Uoxn * Coxn;

% PMOS
Lp   = 0.05;
Vtp  = -0.9056;
Uoxp = 180.2;
Toxp = 2.12e-8;
Coxp = Eox/Toxp;
Kp   = Uoxp * Coxp;

% Step 1
L = 1 * 1e-6;

% Step 2
Cc = 0.22*CL;

% Step 3
I5 = SR*Cc;

% Step 4
K3 = 
Vin = 
Vto3 = 
Vt1 = 
S3 = I5 / ((Kp) * (Vdd - Vin - abs(Vtn) + Vtp)^2);
if S3 < 1
    S3 = 1;
end

% Step 5

