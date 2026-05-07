% =========================================================================
% BINARY
% One time step of the stochastic agent-based model
% Called by sol_MC.m
% "Derivation and quasi-invariant asymptotics of
%  phenotype-structured integro-differential models"
% E. Bernardi, T. Lorenzi, A. Tosin
% =========================================================================

function Vnew = binary(V, Th, MU, Deltat, R, vm, alpha, beta, epsilon, delta, V0)

[N, ~] = size(V);
Vbar   = V;

%% Bernoulli random variables for phenotype changes and compartment switching
Csi1 = binornd(1, MU(1:N/2));
Csi2 = binornd(1, MU(N/2+1:N));

Th1 = binornd(1, Th(1:N/2));
Th2 = binornd(1, Th(N/2+1:N));

%% Random pairing
Ind1 = randperm(N, N/2);
Ind2 = setdiff(1:N, Ind1);

Part1 = Vbar(Ind1, :);
Part2 = Vbar(Ind2, :);

%% Mollifier psi and net proliferation rate r(v)
f_zeta = @(x) 0.5 * (1 + tanh(2 * x ./ (1 - x.^2))) .* (abs(x) < 1);
psi    = @(v) 1 - f_zeta((2*(v - R)/delta) - 1);

Psi1 = zeros(N/2, 1);
Psi2 = zeros(N/2, 1);
for i = 1:N/2
    v1 = Part1(i, 1);
    v2 = Part2(i, 1);
    if v1 > -R - delta && v1 < -R
        Psi1(i) = psi(-v1);
    elseif v1 >= -R && v1 <= R
        Psi1(i) = 1;
    elseif v1 > R && v1 < R + delta
        Psi1(i) = psi(v1);
    end
    if v2 > -R - delta && v2 < -R
        Psi2(i) = psi(-v2);
    elseif v2 >= -R && v2 <= R
        Psi2(i) = 1;
    elseif v2 > R && v2 < R + delta
        Psi2(i) = psi(v2);
    end
end

r1 = (1 - (Part1(:,1) - vm).^2) .* Psi1;
r2 = (1 - (Part2(:,1) - vm).^2) .* Psi2;

%% Spontaneous compartment switching probabilities (death)
D1 = ((1 - r1) .* Deltat) ./ (1 + (1 - r1) .* Deltat);
D2 = ((1 - r2) .* Deltat) ./ (1 + (1 - r2) .* Deltat);

Psi1_bern = binornd(1, D1);
Psi2_bern = binornd(1, D2);

%% Interaction-driven compartment switching (proliferation)
I1  = ones(N/2, 1);
I2  = ones(N/2, 1);
Vp1 = zeros(N/2, 1);
Vp2 = zeros(N/2, 1);

% Pair (0,0)
kk = find(Part1(:,2) + Part2(:,2) == 0);
I1(kk) = 0; I2(kk) = 0;
Vp1(kk) = Part1(kk,1); Vp2(kk) = Part2(kk,1);

% Pair (1,1)
kk = find(Part1(:,2) + Part2(:,2) == 2);
I1(kk) = 1; I2(kk) = 1;
Vp1(kk) = Part1(kk,1); Vp2(kk) = Part2(kk,1);

% Pair (0,1)
kk = find(Part1(:,2) - Part2(:,2) < 0);
I1(kk) = 1; I2(kk) = 1;
Vp1(kk) = Part2(kk,1); Vp2(kk) = Part2(kk,1);

% Pair (1,0)
kk = find(Part2(:,2) - Part1(:,2) < 0);
I1(kk) = 1; I2(kk) = 1;
Vp1(kk) = Part1(kk,1); Vp2(kk) = Part1(kk,1);

%% Structuring-variable switching (phenotype changes)
% Gaussian mutation kernel: v = w + alpha*epsilon^2 + epsilon*sqrt(beta)*Z
VP1 = normrnd(Part1(:,1) + alpha*epsilon^2, epsilon*sqrt(beta));
VP2 = normrnd(Part2(:,1) + alpha*epsilon^2, epsilon*sqrt(beta));

% Individuals in compartment i=0 do not change phenotype
kk = find(Part1(:,2) == 0); VP1(kk) = Part1(kk,1);
kk = find(Part2(:,2) == 0); VP2(kk) = Part2(kk,1);

%% Update compartment and phenotype
Part1(:,2) = (1-Th1).*(1-Psi1_bern).*Part1(:,2) + Th1.*(1-Psi1_bern).*I1 ...
           + Psi1_bern.*(1-Th1).*zeros(N/2,1)    + Psi1_bern.*Th1.*Part1(:,2);

Part1(:,1) = (1-Th1).*(1-Psi1_bern).*((1-Csi1).*Part1(:,1) + Csi1.*VP1) ...
           + Th1.*(1-Psi1_bern).*Vp1 ...
           + Psi1_bern.*(1-Th1).*(V0*ones(N/2,1)) ...
           + Psi1_bern.*Th1.*Part1(:,1);

Part2(:,2) = (1-Th2).*(1-Psi2_bern).*Part2(:,2) + Th2.*(1-Psi2_bern).*I2 ...
           + Psi2_bern.*(1-Th2).*zeros(N/2,1)    + Psi2_bern.*Th2.*Part2(:,2);

Part2(:,1) = (1-Th2).*(1-Psi2_bern).*((1-Csi2).*Part2(:,1) + Csi2.*VP2) ...
           + Th2.*(1-Psi2_bern).*Vp2 ...
           + Psi2_bern.*(1-Th2).*(V0*ones(N/2,1)) ...
           + Psi2_bern.*Th2.*Part2(:,1);

Vbar(Ind1,:) = Part1;
Vbar(Ind2,:) = Part2;
Vnew = Vbar;

end