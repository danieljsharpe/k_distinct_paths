% generic script for a line plot from a .dat file in Matlab
clc, clf;
% absolute path to data file
abs_data_path = '/home/djs244/k_distinct_paths/dummy_fortran/resultscum';
fname = "ngt_cum.dat";
ffmt="%d %f %f %f %f";
col1_idx = 1; % column no. containing x data
col2_idx = 2; % column no. containing y data

% set dimensions of figure and of page
fig_height=400; fig_width=550; % in px
page_height=10; page_width=12.5; % in cm

fid=fopen(fullfile(abs_data_path,fname));
celldata=textscan(fid,ffmt); % read data into cell array
fclose(fid);
xdatacell=celldata(:,col1_idx);
ydatacell=celldata(:,col2_idx);
xdata=cell2mat(xdatacell);
ydata=cell2mat(ydatacell);
%hf = figure(); % octave - what does this do?
plot(xdata,ydata,"b","Linewidth",3.);
hold on;
xlabel('Path number','Interpreter','latex','FontSize',14)
ylabel('Accumulated $k_{AB}^\mathrm{NGT}$','Interpreter','latex','FontSize',14)
%set the x limits and tick intervals here
xlimits=[0,150];
ylimits=[0.,0.2E-02];
xlim([xlimits(1) xlimits(2)])
ylim([ylimits(1) ylimits(2)])
xtickintvl=25;
ytickintvl=0.5E-03;
xticks(xlimits(1):xtickintvl:xlimits(2));
yticks(ylimits(1):ytickintvl:ylimits(2));
ax = gca;
set(ax,'TickDir','out');
set(ax,'TickLabelInterpreter','latex','FontSize',12)
set(gcf,'position',[0,0,fig_width,fig_height]);
set(gcf,'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',[0,0,page_width,page_height]);
% matlab saving syntax
saveas(gcf,fullfile(abs_data_path,'3h_cum_kngt'),'epsc')
saveas(gcf,fullfile(abs_data_path,'3h_cum_kngt'),'png')
% octave saving syntax -- doesn't work!???????????
%print(hf,'threeholepot_pathcosts.pdf','-dpdflatexstandalone');
%system('pdflatex threeholepot_pathcosts');
%open threeholepot_pathcosts.pdf