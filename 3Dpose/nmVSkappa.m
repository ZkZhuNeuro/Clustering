x = linspace(0, 7 / 4 * pi, 8);
DC = rand(50, 1) * 10;
G = rand(50, 1) * 10;
Miu = rand(50, 1) * pi / 4;
Kappa = rand(1, 100) * 18;
figure();
k = [];
n = [];
for i = 1:20
    color = [rand(1), rand(1), rand(1)];
    for j = 1:25
        vec = MakeVonMises4Tuning(linspace(0, 315, 8), [DC(1), G(i), Kappa(j), Miu(i) * 180 / pi]);
        sumvec = sum(vec .* exp(1i * x));
        nm = norm(sumvec) / sum(vec - min(vec));
%         scatter(Kappa, nm, 36, color); hold on
        k = [k ,Kappa(j)];
        n = [n, nm];
    end
end

% xlabel('Kappa');
% % ylabel('Miu');
% ylabel('nm');
scatter(k', n');
[r, p] = corr([k'], [n'], 'Type', 'Spearman');