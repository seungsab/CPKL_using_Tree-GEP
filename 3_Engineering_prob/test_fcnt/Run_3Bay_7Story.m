function y_max = Run_3Bay_7Story(xx)
%% VARIABLE
% xx(1): Width #1 (m) // Norm
% xx(2): Width #2 (m) // Norm

% xx(3): Young's modulus of Floor, Conc (MPa) // Normal
% xx(4): Young's modulus of Steel, Steel (MPa) // Normal

% xx(5): Moment of intertia, floor #1 (m^4) // Normal
% xx(6): Moment of intertia, floor #2 (m^4) // Normal
% xx(7): Moment of intertia, floor #3 (m^4) // Normal
% xx(8): Moment of intertia, floor #4 (m^4) // Normal
% xx(9): Moment of intertia, floor #5 (m^4) // Normal
% xx(10): Moment of intertia, column #1 (m^4) // Normal
% xx(11): Moment of intertia, column #2 (m^4) // Normal
% xx(12): Moment of intertia, column #3 (m^4) // Normal
% xx(13): Moment of intertia, column #4 (m^4) // Normal

% xx(14): Horizontal Load at floor #1 (KN)  // Rayleigh
% xx(15): Horizontal Load at floor #2 (KN)  // Rayleigh
% xx(16): Horizontal Load at floor #3,4,5 (KN)  // Rayleigh

%% Material properties
L1 = 5; L2 = 3.7; % Length of columns (m) 

B(1) = xx(1); % Width #1 (m)
B(2) = xx(2); % Width #2 (m)

% Material properties: Young's modulus
E(1) = xx(3); % Floor // Conc (MPa)
E(2) = xx(4); % Column // Steel (MPa)

% Sectional properties
A1 = 2e-2;  % Floor // Conc // Area (m^2)
A2 = 1e-2;  % Column // Steel // Area (m^2)

% Moment of intertia
% I = 0.03*ones(1,9);
I(1) = xx(5); % floor #1 (m^4)
I(2) = xx(6); % floor #2 (m^4)
I(3) = xx(7); % floor #3 (m^4)
I(4) = xx(8); % floor #4 (m^4)
I(5) = xx(9); % floor #5 (m^4)
I(6) = xx(10); % column #1 (m^4)
I(7) = xx(11); % column #2 (m^4)
I(8) = xx(12); % column #3 (m^4)
I(9) = xx(13); % column #4 (m^4)

% Applied Load
% P = 100*ones(1,5); % Load at floor #1 (KN)
P(1) = xx(14); % Load at floor #1 (KN)
P(2) = xx(15); % Load at floor #2 (KN)
P(3) = xx(16); % Load at floor #3,4,5 (KN)  // Rayleigh


%% STEP 1 - Discretizing the Domain
% [Element Node#1 Node#2 E A I L ANGLE]
Conn = [
    1 1 5 E(2) A2 I(6) L1 90
    2 2 6 E(2) A2 I(7) L1 90
    3 3 7 E(2) A2 I(8) L1 90
    4 4 8 E(2) A2 I(9) L1 90
    5 5 9 E(2) A2 I(6) L2 90
    6 6 10 E(2) A2 I(7) L2 90
    7 7 11 E(2) A2 I(8) L2 90
    8 8 12 E(2) A2 I(9) L2 90
    9 9 13 E(2) A2 I(6) L2 90
    10 10 14 E(2) A2 I(7) L2 90
    11 11 15 E(2) A2 I(8) L2 90
    12 12 16 E(2) A2 I(9) L2 90
    13 13 17 E(2) A2 I(6) L2 90
    14 14 18 E(2) A2 I(7) L2 90
    15 15 19 E(2) A2 I(8) L2 90
    16 16 20 E(2) A2 I(9) L2 90
    17 17 21 E(2) A2 I(6) L2 90
    18 18 22 E(2) A2 I(7) L2 90
    19 19 23 E(2) A2 I(8) L2 90
    20 20 24 E(2) A2 I(9) L2 90
    21 5 6  E(1) A1 I(1) B(1) 0
    22 6 7 E(1) A1 I(1) B(2) 0
    23 7 8 E(1) A1 I(1) B(1) 0
    24 9 10 E(1) A1 I(2) B(1) 0
    25 10 11 E(1) A1 I(2) B(2) 0
    26 11 12 E(1) A1 I(2) B(1) 0
    27 13 14 E(1) A1 I(3) B(1) 0
    28 14 15 E(1) A1 I(3) B(2) 0
    29 15 16 E(1) A1 I(3) B(1) 0
    30 17 18 E(1) A1 I(4) B(1) 0
    31 18 19 E(1) A1 I(4) B(2) 0
    32 19 20 E(1) A1 I(4) B(1) 0
    33 21 22 E(1) A1 I(5) B(1) 0
    34 22 23 E(1) A1 I(5) B(2) 0
    35 23 24 E(1) A1 I(5) B(1) 0
    ];
    
%% STEP 2 - Element Stiffness Matrices
k = {};
for i=1:size(Conn,1)
    k{i}=PlaneFrameElementStiffness(Conn(i,4),Conn(i,5),Conn(i,6),Conn(i,7),Conn(i,8));
end

%% STEP 3 - Assembling the Global Stiffness Matrix
Kg=zeros(24*3,24*3);

for i=1:size(Conn,1)
    Kg=PlaneFrameAssemble(Kg,k{Conn(i,1)},Conn(i,2),Conn(i,3));
end

%% Step 4 ? Applying the Boundary Conditions:
K=Kg(13:end,13:end);

f=zeros(size(K,1),1);
f_IND = 13:12:21*3;

for i=1:size(f_IND,2)
    if i<3
        f(f_IND(i)-12) = P(i);
    else
        f(f_IND(i)-12) = P(end);
    end
    
end

%% 
u=K\f;

y_max = u(58);
end
