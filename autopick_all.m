clear;
close all
wdir='/space/mibhk/hunter/autopick/';
ddir='/space/mibhk/hunter/autopick/data_new';
addpath(wdir); addpath (ddir)

cd (ddir)
tic
unix('ls -d 2*/ > dlist');

%dlist=readlines('dlist');       % readline is not avaiable in 2019 version
%nonEmpty = strlength(dlist) > 0;
%dlist = dlist(nonEmpty);

fid = fopen('dlist','r');
dlist = textscan(fid, '%s');
fclose(fid);
dlist = dlist{:};
% event directory loop
c=1; n = length(dlist);
% event directory loop
for j=1:length(dlist)
    edir=dlist(j);
    cd (char(edir))
    disp (['Doing Event : ', char(edir), ' ==> ', num2str(c), '/', ...
        num2str(n), '...'])
    % station loop
    unix('ls *DPZ*SAC > slist')
    %slist=readlines('slist');
    %  %remove empty string entries
    %nonEmpty = strlength(slist) > 0;
    %slist = slist(nonEmpty);
    fid=fopen('slist', 'r');
    slist = textscan(fid, '%s'); 
    fclose(fid);
    slist = slist{:};
   

    cc=1;nn=length(slist);
    fid = fopen('auto_picked.txt','w');
    for i=1:length(slist)
        close all
        disp (['Event :: ', char(edir), ' Progress ::  ', num2str(cc), ...
            '/', num2str(nn), ' ...'])
        sfile=char(slist(i));
        % get the manual picked arrivals
        waveform=rdsac(sfile);
        if isfield(waveform, 'T1')
            tp_man = waveform.HEADER.T1;
            ts_man = waveform.HEADER.T2;
        else
            tp_man = NaN;
            ts_man = NaN;
        end
        
        % run the autopick to get P and S arrival time picks
        [tp1, ts1,tp2, ts2, t_Pn, t_Sn,coment] = auto_pick_group(sfile);
        TP(i,:) = [tp1, ts1, tp_man, ts_man, t_Pn, t_Sn];
        TP2(i,:) ={tp1, ts1, tp_man, ts_man, t_Pn, t_Sn, coment};
        fprintf(fid, '%s  %7.4f  %7.4f  %7.4f  %7.4f %s\n', string(slist(i)), tp1, ts1, ts2, tp2, string(coment));
    
        cc = cc +1;
    end
    fclose(fid);
    save('picks.mat','TP')
    save('picks_2.mat','TP2')
    cd (ddir)
    c = c+1;
    clear TP
    clear TP2
end

disp('All Finished')

toc

%% plot picks values
%load picks.mat
%load all  picks
wdir ='/space/mibhk/hunter/autopick/data_new';
addpath(wdir)
cd (wdir)
close all
% p1 = load('2023_132_08_33_21\DONE\picks_2.mat');
% p2 = load('2020_113_18_46_42\DONE\picks_2.mat');
% p3 = load('2023_156_16_40_47\DONE\picks_2.mat');
% p4 = load('2023_132_21_58_02\DONE\picks_2.mat'); % event association problem 
% p5 = load('2023_132_00_54_10\DONE\picks_2.mat'); % event superposition
% p6 = load('2023_123_17_16_41\DONE\picks_2.mat'); 
% p7 = load('2023_123_06_20_54\DONE\picks_2.mat');
% p8 = load('2020_111_05_19_11\DONE\picks_2.mat');
% p9 = load('2019_314_22_07_59\DONE\picks_2.mat');
% p10= load('2019_324_06_42_23\DONE\picks_2.mat');

% TP = [p1.TP2; p2.TP2; p3.TP2; p4.TP2; p5.TP2; p6.TP2; p7.TP2; p8.TP2; p9.TP2; p10.TP2];
%TP = [p1.TP; p2.TP; p3.TP; p4.TP; p6.TP; p7.TP; p8.TP; p9.TP; p10.TP];

% load all picks

%TP=[p9.TP2; p10.TP2];
%cell2mat(p1.TP2(:,1));
%TP = [p1.TP; p2.TP; p3.TP; p6.TP];

% filter out the bad event assoication picks
rowsToRemove = any(strcmp('P/S', TP), 2);
% Remove rows with 'ok'
TP(rowsToRemove, :) = [];


tp_dif2 = cell2mat(TP(:,3)) - cell2mat(TP(:,1));
ts_dif2 = cell2mat(TP(:,4)) - cell2mat(TP(:,2));

% remove the NaN values
tp_dif =rmmissing(tp_dif2);
ts_dif = rmmissing(ts_dif2);

tp_good = tp_dif <= 30;
tp_dif = tp_dif(tp_good);
ts_good = ts_dif <= 30;
ts_dif = ts_dif(ts_good);


% figure(4)
% % subplot(2,1,1)
% histogram(tp_dif,40)
% hold on 



%pd = fitdist(tp_dif, 'beta', 'lower', -Inf, 'upper', Inf);
% %Plot the fitted beta distribution curve




fig = figure('position',[500, 300, 800, 600]);

pd = fitdist(tp_dif, "Normal");
%x = linspace(min(tp_dif), max(tp_dif), 100);
x = linspace(min(tp_dif), -min(tp_dif), 100);
y = pdf(pd, x);
%h = histogram(tp_dif,40,'FaceColor',[0 0.4470 0.7410]);
h = histogram(tp_dif,50,'Normalization','pdf','FaceColor', ...
    [0 0.4470 0.7410], 'FaceAlpha',0.5);

hold on
% Add a second y-axis on the right
%yyaxis right;
plot(x, y, 'LineWidth', 2,'LineStyle','-','Color',[0.6350 0.0780 0.1840]);
hold on
line([pd.mean+pd.sigma, pd.mean+pd.sigma ], [0, pdf(pd,(pd.mean+pd.sigma))], ...
    'LineWidth',1.5,'Color','g', 'LineStyle', '--')
hold on
line([pd.mean-pd.sigma, pd.mean-pd.sigma ], [0, pdf(pd,(pd.mean-pd.sigma))], ...
    'LineWidth',1.5,'Color','g', 'LineStyle', '--')
hold on
line([pd.mean, pd.mean ], [0, pdf(pd,pd.mean)], ...
    'LineWidth',1.5,'Color','k', 'LineStyle', '--')
hold on
str1 = ['Mean: ', sprintf('%04.3f',pd.mean),' s'];
str2 = ['Std: ', sprintf('%04.3f',pd.sigma),' s'];
str3 = ['N : ', sprintf('%d', length(pd.InputData.data))];
annotation('textbox', [0.7, 0.82, 0.1, 0.1], ...
    'String',str1, 'FontSize', 12, 'Color', 'k', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'EdgeColor','None');
hold on

annotation('textbox', [0.7, 0.78, 0.1, 0.1], ...
    'String',str2, 'FontSize', 12, 'Color', 'k', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'EdgeColor','None');

hold on
annotation('textbox', [0.7, 0.74, 0.1, 0.1], ...
    'String',str3, 'FontSize', 12, 'Color', 'k', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'EdgeColor','None');


ylabel('Probabilty Distribution','FontSize',14,'FontWeight','bold')
xlabel('Pick Misfit (s)','FontSize',14,'FontWeight','bold')
title('Manual vs Autopick Misfit Distribution [P wave]', 'FontSize',16, ...
    'FontWeight','bold')

filename='P_misfit.png';
exportgraphics(gcf,filename,'Resolution',300)

%% %% plot S wave misfit

fig2 = figure('position',[500, 300, 800, 600]);

pd = fitdist(ts_dif,'Normal');
%x = linspace(min(ts_dif), max(ts_dif), 100);
x = linspace(min(ts_dif), -min(ts_dif), 100);
y = pdf(pd, x);
%h = histogram(tp_dif,40,'FaceColor',[0 0.4470 0.7410]);
h = histogram(ts_dif,50,'Normalization','pdf','FaceColor', ...
    [0 0.4470 0.7410], 'FaceAlpha',0.5);

hold on
% Add a second y-axis on the right
%yyaxis right;
plot(x, y, 'LineWidth', 2,'LineStyle','-','Color',[0.6350 0.0780 0.1840]);
hold on
line([pd.mean+pd.sigma, pd.mean+pd.sigma ], [0, pdf(pd,(pd.mean+pd.sigma))], ...
    'LineWidth',1.5,'Color','g', 'LineStyle', '--')
hold on
line([pd.mean-pd.sigma, pd.mean-pd.sigma ], [0, pdf(pd,(pd.mean-pd.sigma))], ...
    'LineWidth',1.5,'Color','g', 'LineStyle', '--')
hold on
line([pd.mean, pd.mean ], [0, pdf(pd,pd.mean)], ...
    'LineWidth',1.5,'Color','k', 'LineStyle', '--')
hold on
str1 = ['Mean: ', sprintf('%04.3f',pd.mean),' s'];
str2 = ['Std: ', sprintf('%04.3f',pd.sigma),' s'];
str3 = ['N : ', sprintf('%d', length(pd.InputData.data))];
annotation('textbox', [0.7, 0.82, 0.1, 0.1], ...
    'String',str1, 'FontSize', 12, 'Color', 'k', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'EdgeColor','None');
hold on

annotation('textbox', [0.7, 0.78, 0.1, 0.1], ...
    'String',str2, 'FontSize', 12, 'Color', 'k', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'EdgeColor','None');

hold on
annotation('textbox', [0.7, 0.74, 0.1, 0.1], ...
    'String',str3, 'FontSize', 12, 'Color', 'k', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'EdgeColor','None');


ylabel('Probabilty Distribution','FontSize',14,'FontWeight','bold')
xlabel('Pick Misfit (s)','FontSize',14,'FontWeight','bold')
title('Manual vs Autopick Misfit Distribution [S wave]', 'FontSize',16, ...
    'FontWeight','bold')

filename='S_misfit.png';
exportgraphics(gcf,filename,'Resolution',300)
% figure(6)
% % beta distribution
% 
% tpdif_nor = (tp_dif - min(tp_dif)) / (max(tp_dif) - min(tp_dif));
% tpdif_back = tpdif_nor * (max(tp_dif) - min(tp_dif)) + min(tp_dif);
% %histogram(tpdif_nor,40)
% h = histogram(tp_dif,50,'Normalization','pdf','FaceColor',[0 0.4470 0.7410]);
% hold on
% pd2 = fitdist(tpdif_nor, 'Beta');
% x2 = linspace(min(tpdif_nor), max(tpdif_nor), 100);
% y2 = pdf(pd2, x2);
% plot(x2, y2,'r', 'LineWidth', 2);


% 
% % histfit(tpdif_nor,40,'beta');
% title('Manual vs Autopick [P-wave]','FontSize',14,'FontWeight','bold')
% subplot(2,1,2)
% histogram(ts_dif,40)
% hold on
% h2 = histfit(ts_dif,40,'beta');
% title('Manual vs Autopick [S-wave]','FontSize',14, 'FontWeight', 'bold')
% filename='man_vs_auto.png';
% exportgraphics(gcf,filename,'Resolution',300)

