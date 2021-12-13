clear
close all

%% Load data
rgb_shifts = struct();

%%%%%%%%%%%%%%%%%%%%%%%%%%%CHANGE TO YOUR PATH%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rgb_shifts.analysis_path = ''; %Main data path: '/Example_Input/'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outdir = [rgb_shifts.analysis_path filesep 'RGB PSF/'];
mkdir(outdir)
cd(outdir)

%%%RECONSTRUCTIONS ARE INPUTS HERE (THIS STRUCT CONTAINS INFO ABOUT THE PSF)
red_psf_loc = [rgb_shifts.analysis_path filesep 'red' filesep 'recon_out' filesep];
cd(red_psf_loc);
psf = dir('*.mat');
red_psf = psf.name;

green_psf_loc = [rgb_shifts.analysis_path filesep 'green' filesep 'recon_out' filesep];
cd(green_psf_loc);
psf = dir('*.mat');
green_psf = psf.name;

blue_psf_loc = [rgb_shifts.analysis_path filesep 'blue' filesep 'recon_out' filesep];
cd(blue_psf_loc);
psf = dir('*.mat');
blue_psf = psf.name;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([rgb_shifts.analysis_path filesep 'red' filesep 'recon_out' filesep red_psf]);
rgb_shifts.red_psf = image_params;
load([rgb_shifts.analysis_path filesep 'green' filesep 'recon_out' filesep green_psf]);
rgb_shifts.green_psf = image_params;
load([rgb_shifts.analysis_path filesep 'blue' filesep 'recon_out' filesep blue_psf]);
rgb_shifts.blue_psf = image_params;

step_size = 0.050;

%% Populate shifts structure

[rgb_shifts.red_total_shifts,rgb_shifts.red_off_ax,rgb_shifts.red_off_ax_av,rgb_shifts.red_row_shifts,rgb_shifts.red_col_shifts] = ...
    PlotShiftsOverDistance(rgb_shifts.red_psf,step_size);

[rgb_shifts.blue_total_shifts,rgb_shifts.blue_off_ax,rgb_shifts.blue_off_ax_av,rgb_shifts.blue_row_shifts,rgb_shifts.blue_col_shifts] = ...
    PlotShiftsOverDistance(rgb_shifts.blue_psf,step_size);

[rgb_shifts.green_total_shifts,rgb_shifts.green_off_ax,rgb_shifts.green_off_ax_av,rgb_shifts.green_row_shifts,rgb_shifts.green_col_shifts] = ...
    PlotShiftsOverDistance(rgb_shifts.green_psf,step_size);

%% Initialize Parameters
n_lenses = size(rgb_shifts.green_total_shifts,1);
n_off_ax = n_lenses - 1;
total_shifts = size(rgb_shifts.green_total_shifts,2);

%%%% THIS WILL DEPEND ON IF THE PSF WAS TAKEN CLOSE TO FAR OR FAR TO CLOSE
% close to far
%xax = 0:step_size:(step_size*total_shifts-step_size);
% far to close
xax = (step_size*total_shifts-step_size):-step_size:0;
%% RGB shift per lens
figure
ha = tight_subplot(2,3,0.1,0.05,0.05);
p = 1;
for i=1:n_off_ax
    set(gcf,'CurrentAxes',ha(i))
    plot(xax,rgb_shifts.red_off_ax(i,:),'-r');
    hold on
    plot(xax, rgb_shifts.green_off_ax(i,:),'-g');
    plot(xax, rgb_shifts.blue_off_ax(i,:),'-b');
    box on
    set(gca,'fontsize',35) % 10
    set(gca,'linewidth',3) % 3
    set(findobj(gca,'type','line'),'linew',2)
    hold off
    p = p+1;
end
set(gcf,'CurrentAxes', ha(3))
legend({'Red','Green','Blue'},'Fontsize',25)

%%% RGB Row and Col shifts %%%

%% Row shifts
figure
hold on
for i=1:n_lenses
    plot(xax,rgb_shifts.red_row_shifts(i,:),'.r')
    plot(xax,rgb_shifts.green_row_shifts(i,:),'.g')
    plot(xax,rgb_shifts.blue_row_shifts(i,:),'.b')
end  
%legend({'1', '2', '3', '4', '5', '6', '7'},'Orientation','horizontal','Location','Best');
box on
set(gca,'fontsize',35) % 10
set(gca,'linewidth',3) % 3
set(findobj(gca,'type','line'),'linew',3)

%% Column shifts
figure
hold on
for i=1:n_lenses
    plot(xax,rgb_shifts.red_col_shifts(i,:),'.r')
    plot(xax,rgb_shifts.green_col_shifts(i,:),'.g')
    plot(xax,rgb_shifts.blue_col_shifts(i,:),'.b')
end  
%legend({'1', '2', '3', '4', '5', '6', '7'},'Orientation','horizontal','Location','Best');
box on
set(gca,'fontsize',35) % 10
set(gca,'linewidth',2) % 3
set(findobj(gca,'type','line'),'linew',2)

%% Column Shifts 2
figure
hold on
for i=1:n_lenses
    plot(xax,rgb_shifts.red_col_shifts(i,:),'.r','MarkerSize',10)
    plot(xax,rgb_shifts.green_col_shifts(i,:),'.g','MarkerSize',10)
    plot(xax,rgb_shifts.blue_col_shifts(i,:),'.b','MarkerSize',10)
end  
%legend({'1', '2', '3', '4', '5', '6', '7'},'Orientation','horizontal','Location','Best');
box on
set(gca,'fontsize',35) % 10
set(gca,'linewidth',5) % 3
set(findobj(gca,'type','line'),'linew',10)

%% RGB total shift

figure
plot(xax, rgb_shifts.red_off_ax_av,'-r')
hold on
plot(xax, rgb_shifts.green_off_ax_av,'-g')
plot(xax, rgb_shifts.blue_off_ax_av, '-b')
box on
set(gca,'fontsize',35) % 10
set(gca,'linewidth',3) % 3
set(findobj(gca,'type','line'),'linew',2)
legend({'Red','Green','Blue'},'Fontsize',25)

%% Magnification Calculation
pitch = rgb_shifts.green_psf.lens_pitch;
pixel_D = rgb_shifts.green_psf.Pixel_D;

% Red
rgb_shifts.red_mag = Mag(rgb_shifts,'red',step_size,15,3,pitch,pixel_D);
% Green
rgb_shifts.green_mag = Mag(rgb_shifts,'green',step_size,15,3,pitch,pixel_D);
% Blue
rgb_shifts.blue_mag = Mag(rgb_shifts,'blue',step_size,15,3,pitch,pixel_D);

%% Fit to system magnification function to get lens spacing inside endoscope
zrange = xax;
%zrange = step_size:step_size:size(rgb_shifts.green_total_shifts,2)*step_size;

figure('Position',[360,82,560,616])
stdshade(rgb_shifts.green_mag.all_lens_m,0.3,'green',zrange);
hold on
p_g = plot(zrange,rgb_shifts.green_mag.avg_m,'linewidth',1.5,'color','green')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Magnification Modeling

% All lenses all colors magnification plot, average, curvefit
%zrange = step_size:step_size:size(rgb_shifts.red_total_shifts,2)*step_size;
zrange = xax;

figure('Position',[360,82,560,616])
stdshade(rgb_shifts.red_mag.all_lens_m,0.3,'red',zrange);
hold on
stdshade(rgb_shifts.blue_mag.all_lens_m,0.3,'blue',zrange);
stdshade(rgb_shifts.green_mag.all_lens_m,0.3,'green',zrange);

p_r = plot(zrange,rgb_shifts.red_mag.avg_m,'linewidth',1.5,'color','red')

p_b = plot(zrange,rgb_shifts.blue_mag.avg_m,'linewidth',1.5,'color','blue')

p_g = plot(zrange,rgb_shifts.green_mag.avg_m,'linewidth',1.5,'color','green')


rgb_shifts.avg_m_allcolor = mean([rgb_shifts.red_mag.avg_m;rgb_shifts.green_mag.avg_m;rgb_shifts.blue_mag.avg_m]);
rgb_shifts.avg_m_allcorlor_f = fitmag(zrange',rgb_shifts.avg_m_allcolor',16,5)
pf = plot(zrange,rgb_shifts.avg_m_allcorlor_f(zrange),'--','linewidth',2.5,'color','black')

rgb_shifts.avg_m_allcolor_fit = rgb_shifts.avg_m_allcorlor_f(zrange);

title(['Lens spacing (a) = ' num2str(rgb_shifts.avg_m_allcorlor_f.a) newline 'PSF Offset (o) = ' num2str(rgb_shifts.avg_m_allcorlor_f.o)])

box on
set(gca,'fontsize',35) % 10
set(gca,'linewidth',3) % 3
legend([p_r,p_g,p_b,pf],'M_R','M_G','M_B','fit')

%% Average RGB fits with SET lens spacing

% do the fits
%y = magfit_setspacing(zrange,rgb_shifts.avg_m_allcorlor_f.o)
rgb_shifts.red_m_fixed_ls_f = fitmag_setspacing(zrange', rgb_shifts.red_mag.avg_m',rgb_shifts.avg_m_allcorlor_f.o)
rgb_shifts.green_m_fixed_ls_f = fitmag_setspacing(zrange', rgb_shifts.green_mag.avg_m',rgb_shifts.avg_m_allcorlor_f.o)
rgb_shifts.blue_m_fixed_ls_f = fitmag_setspacing(zrange', rgb_shifts.blue_mag.avg_m',rgb_shifts.avg_m_allcorlor_f.o)

% red
%figure('Position',[360,57,728,641])
figure('Position',[360,82,560,616])
stdshade(rgb_shifts.red_mag.all_lens_m,0.3,'red',zrange);
hold on
p_r2 = plot(zrange,rgb_shifts.red_mag.avg_m,'linewidth',1.5,'color','red')
p_rf = plot(zrange,rgb_shifts.red_m_fixed_ls_f(zrange),'--black','linewidth',2)

rgb_shifts.red_m_fixed_ls_fit = rgb_shifts.red_m_fixed_ls_f(zrange);

box on
title(['Lens spacing (a) = ' num2str(rgb_shifts.avg_m_allcorlor_f.a) newline 'PSF Offset (o) = ' num2str(rgb_shifts.red_m_fixed_ls_f.o)])

set(gca,'fontsize',35) % 10
set(gca,'linewidth',2) % 3
%set(p1,'LineWidth',4)
%set(p2,'LineWidth',2)
%set(findobj(gca,'type','line'),'linew',2)
xlabel('')
ylabel('');
ylim([0 0.37])
legend([p_r2,p_rf],'Experimental', 'Fit','Fontsize',25)

% green
figure('Position',[360,82,560,616])
stdshade(rgb_shifts.green_mag.all_lens_m,0.3,'green',zrange);
hold on
p_g2 = plot(zrange,rgb_shifts.green_mag.avg_m,'linewidth',1.5,'color','green')
p_gf = plot(zrange,rgb_shifts.green_m_fixed_ls_f(zrange),'--black','linewidth',2)

rgb_shifts.green_m_fixed_ls_fit = rgb_shifts.green_m_fixed_ls_f(zrange);

box on
title(['Lens spacing (a) = ' num2str(rgb_shifts.avg_m_allcorlor_f.a) newline 'PSF Offset (o) = ' num2str(rgb_shifts.green_m_fixed_ls_f.o)])

set(gca,'fontsize',35) % 10
set(gca,'linewidth',2) % 3
%set(p1,'LineWidth',4)
%set(p2,'LineWidth',2)
%set(findobj(gca,'type','line'),'linew',2)
xlabel('')
ylabel('');
ylim([0 0.37])
legend([p_g2,p_gf],'Experimental', 'Fit','Fontsize',25)


% blue
figure('Position',[360,82,560,616])
stdshade(rgb_shifts.blue_mag.all_lens_m,0.3,'blue',zrange);
hold on
p_b2 = plot(zrange,rgb_shifts.blue_mag.avg_m,'linewidth',1.5,'color','blue')
p_bf = plot(zrange,rgb_shifts.blue_m_fixed_ls_f(zrange),'--black','linewidth',2)

rgb_shifts.blue_m_fixed_ls_fit = rgb_shifts.blue_m_fixed_ls_f(zrange);

box on
title(['Lens spacing (a) = ' num2str(rgb_shifts.avg_m_allcorlor_f.a) newline 'PSF Offset (o) = ' num2str(rgb_shifts.blue_m_fixed_ls_f.o)])

set(gca,'fontsize',35) % 10
set(gca,'linewidth',2) % 3
%set(p1,'LineWidth',4)
%set(p2,'LineWidth',2)
%set(findobj(gca,'type','line'),'linew',2)
xlabel('')
ylabel('');
ylim([0 0.37])
legend([p_b2,p_bf],'Experimental', 'Fit','Fontsize',25)


%% SAVE FIGS

FolderName = outdir;   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = num2str(get(FigHandle, 'Number'));
  savefig(FigHandle, [FolderName filesep FigName '.fig'])%fullfile(FolderName, FigName, '.fig'));
  saveas(FigHandle, [FolderName filesep FigName '.png'])
end

save([rgb_shifts.analysis_path filesep 'RGB PSF/rgb_shifts.mat'],'rgb_shifts')


%%
function f = fitmag(x,y,a,o)
    ft = fittype( 'magfit(x,a,o)');
    f = fit(x,y,ft,'StartPoint',[a,o]);
end


function f = fitmag_setspacing(x,y,o)
    ft = fittype( 'magfit_setspacing(x,o)');
    f = fit(x,y,ft,'StartPoint',o);
end

function mag = Mag(rgb_shifts,color,step_size,grin_to_relay,offset,pitch,pixel_D)
    
    fit_params = struct()
    
    image_params = rgb_shifts.([color '_psf'])
    
    switch color
        case 'red'
            linecolor = 'r';
        case 'blue'
            linecolor = 'b';
        case 'green'
            linecolor = 'g';
    end
    
    % Off axis shifts and average
    off_ax = rgb_shifts.([color '_off_ax'])
    size(off_ax)
    off_ax_av = rgb_shifts.([color '_off_ax_av']);      
    zsteps = length(off_ax_av);
    zstepsize = step_size;
    
    % Axial distance
    %%%%%% THIS WILL DEPEND ON IF THE PSF IS TAKEN FAR TO CLOSE OR CLOSE TO
    %%%%%% FAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % close to far
    %zrange = zstepsize:zstepsize:zsteps*zstepsize;
    % far to close
    zrange = zsteps*zstepsize:-zstepsize:zstepsize;
    
%     pitch = image_params.lens_pitch;
    lens_D = image_params.circle_radius*2;
    obj_pix_size = 1/lens_D;
    %im_pix_size = image_params.Pixel_D;
    im_pix_size = pixel_D;
    n_lenses = size(rgb_shifts.green_total_shifts,1);
    n_off_ax = n_lenses - 1;
    
    % Fit Params
    x = zrange';
    ft = fittype( 'magfit(x,a,o)');
    
    %% Plot all lenses and average
    figure
    p1 = plot(x,off_ax,linecolor,'Linewidth',1)
    hold on
    p2 = plot(x,off_ax_av,'black','Linewidth',3)
    l = legend([p1(1),p2],'Individual lenses','Average')
    set(l,'Fontsize',25)
    box on
    set(gca,'fontsize',35) % 10
    set(gca,'linewidth',3) % 3
    hold off
    
    %% Plot curve fit for each lens
    
    % magnification
    m = (off_ax.*im_pix_size)/pitch;
    m1 = m';
    fit_params.all_lens_m = m;
    %total_shifts = size(rgb_shifts.red_total_shifts,2);
    %xax = 0:step_size:(step_size*total_shifts-step_size);
    figure('units','normalized','outerposition',[0 0 1 1],'Name',color);
    ha = tight_subplot(2,3,[0.2 0.1],[0.05,0.15],0.07);
    for a=1:n_off_ax    
        % fit 
        f_i = fit(x,m1(:,a),ft,'StartPoint',[grin_to_relay,offset]);
        fit_params.lensfit{a,1} = f_i;
        
        % plot
        sprintf(num2str(a))
        set(gcf,'CurrentAxes',ha(a))
        p1 = plot(f_i,x,m1(:,a));
        box on
        set(gca,'fontsize',20) % 10
        set(gca,'linewidth',2) % 3
        set(p1(1),'LineWidth',4)
        set(p1(2),'LineWidth',1)
        %set(findobj(gca,'type','line'),'linew',2)
        title(['Lens spacing = ' num2str(f_i.a) newline 'Offset = ' num2str(f_i.o)])
        xlabel('')
        ylabel('');
        if a ~= 3
            legend('off')
        end
    end

    %% Plot curve fit for all lenses

     magnification = (off_ax_av.*im_pix_size)/pitch;
     y = magnification';
     f = fit(x,y,ft,'StartPoint',[15,3]);

     % Plot
    %figure('Position',[360,57,728,641])
    %figure('Position',[229,138,728,391])
    figure('Position',[360,285,444,413])
    stdshade(m,0.2,linecolor,zrange);
    hold on
    p1 = plot(x,y,linecolor,'LineWidth',2);
    p2 = plot(x,f(x),'black--');
    box on
    title(['Lens spacing = ' num2str(f.a) newline ' Offset = ' num2str(f.o)])
    set(gca,'fontsize',35) % 10
    set(gca,'linewidth',2) % 3
    %set(p1,'LineWidth',4)
    set(p2,'LineWidth',2)
    %set(findobj(gca,'type','line'),'linew',2)
    xlabel('')
    ylabel('');
    ylim([0 0.37])
    legend([p1,p2],'Experimental', 'Fit','Fontsize',25)

    % Housekeeping
    fit_params.avg_m = magnification;
    fit_params.avg_fit = y; 
    mag = fit_params;

    %     plot(f,x,y)
    %     graph_beautify(gca)
    %     
    %     fit_params.x = x;
    %     fit_params.y = y;
    %     fit_params.f = f;

    %mag.off_ax_shifts = off_ax;
    %mag.off_ax_av_shifts = off_ax_av;

    %image_params.mag = mag;
    %image_params.fit = fit_params;
end