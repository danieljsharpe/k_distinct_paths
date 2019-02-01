clc, clf;
n_paths = 15; % choose no. of paths to plot (no. of Epath files)
max_path_length = 80;
abs_data_path = '/home/djs244/HG_WC/kdp_py_results'; % absolute path to data (Epath) files
paths = zeros(max_path_length,n_paths,"double");
sp_idx=1:1:max_path_length;
for i = 1:n_paths
    fid = fopen(fullfile(abs_data_path,strcat("Epath.",string(i),".txt")));
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
    plot(sp_idx(ix),paths(ix,i),"LineWidth",1.0); % plot only these corresponding values
    hold on;
end
xlabel('Stationary point','Interpreter','latex','FontSize',14)
ylabel('Potential energy / kcal mol$^{-1}$','Interpreter','latex','FontSize',14)
ax = gca;
set(ax,'TickDir','out');
set(ax,'TickLabelInterpreter','latex','FontSize',12)
saveas(gcf,fullfile(abs_data_path,'hg_wc_15paths'),'epsc')
saveas(gcf,fullfile(abs_data_path,'hg_wc_15paths'),'png')