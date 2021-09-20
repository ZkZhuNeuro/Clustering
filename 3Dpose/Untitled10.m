% DC = 1;
% G = (rand(1, 100)) .^ 2 * 10;
% Kappa = 1;
% Miu = 0;
x = linspace(0, 7 / 4 * pi, 8);
% VM = @(x, DC, G, Kappa, Miu) DC + G .* sum(exp(Kappa .* (cos(x - Miu) - 1)));
% vmSum = arrayfun(@(k) VM(x, DC, G, k, Miu), Kappa);
% vmSum = DC + G .* sum(exp(Kappa .* (cos(x - Miu) - 1)));
% vmPeak = DC + G;
% normMag = vmPeak ./ vmSum;
% figure();
% scatter(Kappa, normMag)
DC = 0;
G = rand(50, 1) * 10;
Miu = rand(50, 1) * 10;
Kappa = rand(1, 100) * 10;


% v = (DC + G) / (8 * DC + G * sum(exp(Kappa * (cos(x - Miu) - 1))));
figure();
for i = 1:100
    %     for j = 1:30
    vec = MakeVonMises4Tuning(linspace(0, 315, 8), [DC * ones(50, 1), ...
        G, Kappa(i) * ones(50, 1), Miu * 180 / pi]);
    nvec = (sqrt(real((sum(vec .* exp(1i * linspace(0, 7/4 * pi, 8)), 2)))) .^2 +...
        (imag((sum(vec .* exp(1i * linspace(0, 7/4 * pi, 8)), 2)))) .^2) ./ sum(vec, 2);
    %     S = @(x, Kappa, Miu) sum(exp(Kappa * (cos(x - Miu) - 1)));
    %     s = arrayfun(@(k) S(x, k, Miu), Kappa);
    %     nm = (DC + G(1, i)) ./ (8 * DC + G(1, i) * s);
    scatter3(Kappa(i) * ones(50, 1), G, nvec); hold on
    
    %     end
    % hold on
    % Kappa = 0.0000001;
    % S = @(x, Kappa, Miu) sum(exp(Kappa * cos(x - Miu) - 1));
    % s2 = arrayfun(@(m) S(x, Kappa, m), Miu);
    % scatter(Miu, s2);
    % hold on
    % ylim([0,5]);
end
xlabel('Kappa');
ylabel('G');
zlabel('nm');
