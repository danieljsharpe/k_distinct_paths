clc, clf;
n_paths = 150; % choose no. of paths to plot (no. of Epath files)
%choose one
path_nos=linspace(1,n_paths,n_paths);
%path_nos=[1 10000 20000 50000]; % set n_paths equal to len(path_nos)

%fig_height=400; fig_width=550; % in px
%page_height=10; page_width=12.5; % in cm

max_path_length = 100;
% absolute path to data (Epath) files
abs_data_path = '/home/djs244/k_distinct_paths/dummy_tests3';
paths = zeros(max_path_length,n_paths,"double");
sp_idx=1:1:max_path_length;
my_cmap=winter(n_paths); % choose colourmap, e.g. winter, Rdbu, brg
for i = 1:n_paths % use this FOR statement if set n_paths    
    fid = fopen(fullfile(abs_data_path,strcat("Epath.",string(path_nos(i)),".txt")));
    % fgetl(fid); % read and discard the first line (i.e. if have header)
    % ignore first and third columns, read second column:
    data = fscanf(fid,"%*d %f %*d");
    % data = readtable(strcat("Epath.",string(i))) % file extension not
    % valid
    fclose(fid);
    d = size(data);
    size(padarray(data,max_path_length-d(1),0.,"post"));
    size(paths(:,i));
    paths(:,i) = padarray(data,max_path_length-d(1),0.,"post");
    % CHEATING - changing the TS energies
    % for j = 1:max_path_length
    %   if  paths(j,i) > 20.0
    %       paths(j,i) = paths(j,i) - 15.0
    %   end
    % end
    % N.B. matlab will plot zeros (but not NaNs). Since we have zeros...
    ix=(paths(:,i)~=0); % get indices of non-zero values only
    plot(sp_idx(ix),paths(ix,i),"Color",my_cmap(i,:),"LineWidth",1.5); % plot only these corresponding values
    hold on;
end
xlabel('Stationary point on pathway','Interpreter','latex','FontSize',14)
%ylabel('Potential energy / kcal mol$^{-1}$','Interpreter','latex','FontSize',14)
ylabel('Potential energy $V$','Interpreter','latex','FontSize',14);
%% use legend if plotting specific data sets
%lgd = legend('1','10000','20000','50000','Interpreter','latex','FontSize',6,'location','northeast');
%legend('boxoff');
%title(lgd,'Path number','Interpreter','latex','FontSize',6);
%% use colourbar if plotting range of data sets
colormap(my_cmap);
caxis([1 n_paths]); % rescale colormap
cbar = colorbar;
cbar.Limits = [1 n_paths];
cbar.Ticks = [1 25 50 75 100 125 150];
cbar.TickDirection = 'out';
cbar.FontSize = 12;
cbar.TickLabelInterpreter = 'latex';
cbar.Label.String = 'Path number';
cbar.Label.FontSize = 14;
cbar.Label.Interpreter = 'latex';
cbar.Location = 'northoutside';

ax = gca;
set(ax,'TickDir','out');
set(ax,'TickLabelInterpreter','latex','FontSize',12)
set(gcf,'Position',[0,0,fig_width,fig_height]);
set(gcf,'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',[0,0,page_width,page_height]);
saveas(gcf,fullfile(abs_data_path,'kdp_3h_epaths'),'epsc')
saveas(gcf,fullfile(abs_data_path,'kdp_3h_epaths'),'png')