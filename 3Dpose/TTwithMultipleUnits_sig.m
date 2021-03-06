%%
close all; clear all; clc
%% Load all the data from selected animals' selected areas. 

Subj = {'Lue'}; % 'Lue', 'Frank', 'Warren'
RecordingArea = 'CIP';
names = fileSeeker(Subj, RecordingArea);

%% Save the filenames of tetrodes with multiple units. 
[tetrodes] = fileFilter_tetrodeWithMultipleUnits(names);

%% Analyze the units data, and do the pariwise comparison. 
[pwc, mag, kappa] = pairWiseComp_twmu_sig(tetrodes);

figure();
scatter(mag, kappa); hold on
% scatter(magNs, kappaNs); 
[r_mk, p_mk] = corr(mag', kappa', 'Type', 'Spearman');
s_r_mk = num2str(r_mk);
s_p_mk = num2str(p_mk);
text(0.1, 17, ['r = ', s_r_mk]);
text(0.1, 16, ['p = ', s_p_mk]);
ylim([0, 18]);
xlim([0, 1]);
xlabel('Normalized Magnitude');
ylabel('Kappa');
title('Normalized magnitudes VS von Mises kappas');

%% Plot the scatter plots. 
figure();
subplot(2, 2, 1);

s1 = scatter(pwc.dir_x, pwc.dir_y); 
s1.MarkerFaceColor = [0 0.5 0.5];hold on
% cosDir = cosd(pwc.dir_x - pwc.dir_y);
% sigmaDir = std(cosDir);
% [hDir, pDir] = ztest(cosDir, 0, sigmaDir);
% s_pDir = num2str(pDir);
% text(30, 500, ['p(sumVector) = ', s_pDir]);

s2 = scatter(pwc.dir_vm_x, pwc.dir_vm_y);
s2.MarkerFaceColor = [0.75 0 0];
% cosDirVm = cosd(pwc.dir_vm_x - pwc.dir_vm_y);
% sigmaDirVm = std(cosDirVm);
% [hDirVm, pDirVm] = ztest(cosDirVm, 0, sigmaDirVm);
% s_pDirVm = num2str(pDirVm);
% text(30, 450, ['p(vonMisesMiu) = ', s_pDirVm]);

plot([0, 360], [0, 360], '-k');
xticks(linspace(0, 360, 13));
yticks(linspace(-180, 540, 25));
xlim([0 360]);
ylim([-180 540]);
title('Differences between preferred directions');

 
subplot(2, 2, 2);
scatter(pwc.sdi_x, pwc.sdi_y); hold on
plot([0, 1], [0, 1], '-k');
[r_sdi, p_sdi] = corr(pwc.sdi_x', pwc.sdi_y');
s_r_sdi = num2str(r_sdi);
s_p_sdi = num2str(p_sdi);
text(0.1, 0.9, ['r = ', s_r_sdi]);
text(0.1, 0.8, ['p = ', s_p_sdi]);
title('Distribution of saccade discrimination indices');
xlim([0 1]);
ylim([0 1]);

subplot(2, 2, 3);
scatter(pwc.mag_x, pwc.mag_y); hold on
plot([0, 1], [0, 1], '-k');
[r_mag, p_mag] = corr(pwc.mag_x', pwc.mag_y');
s_r_mag = num2str(r_mag);
s_p_mag = num2str(p_mag);
text(0.1, 0.9, ['r = ', s_r_mag]);
text(0.1, 0.8, ['p = ', s_p_mag]);
title('Distribution of normalized magnitudes');
xlim([0 1]);
ylim([0 1]);

subplot(2, 2, 4);
scatter(pwc.bw_x, pwc.bw_y); hold on
plot([0, 18], [0, 18], '-k');
[r_bw, p_bw] = corr(pwc.bw_x', pwc.bw_y');
s_r_bw = num2str(r_bw);
s_p_bw = num2str(p_bw);
text(1, 17, ['r = ', s_r_bw]);
text(1, 16, ['p = ', s_p_bw]);
title('Distribution of von Mises kappas');
xlim([0 18]);
ylim([0 18]);

%% Plot the differences between preferred directions. 

figure();
subplot(1, 2, 1);
edges = linspace(0, 180, 7);
d_dir = [];
for i = 1:size(pwc.d_dir, 2)
    if pwc.d_dir(i) > 180
        d_dir(i) = 360 - pwc.d_dir(i);
    else 
        d_dir(i) = pwc.d_dir(i);
    end
    
    if pwc.d_dir_vm(i) > 180
        d_dir_vm(i) = 360 - pwc.d_dir_vm(i);
    else
        d_dir_vm(i) = pwc.d_dir_vm(i);
    end
end
histogram(d_dir, edges); hold on
title('Differences between preferred directions');
xlabel('Difference (degree)');
ylabel('Count');
histogram(d_dir_vm, edges);
xlim([0 180]);
xticks(linspace(0, 180, 7));
legend('Sum vectors', 'Von Mises mius');
% ylim([0 180]);

subplot(1, 2, 2);
% edges_rose = linspace(0, pi, 13);
% disp(pwc.d_dir * 180 / pi)
% disp(pwc.d_dir_vm)
edges_rose = linspace(0, 2 * pi, 13);
polarhistogram(pwc.d_dir * pi / 180,'BinEdges', ...
    edges_rose + 15 / 180 * pi); hold on
polarhistogram(pwc.d_dir_vm * pi / 180,'BinEdges', ...
    edges_rose + 15 / 180 * pi); 
title('Differences between preferred directions');
legend('Sum vectors', 'Von Mises mius','Location','northeast');