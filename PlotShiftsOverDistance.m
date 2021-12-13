function [total_shifts,off_ax,off_ax_av,row_shifts,col_shifts] = PlotShiftsOverDistance(image_params,step_size)
    %step_size = 0.05;
    shifts = image_params.shifts_normalized;  

    %% Compute total diagonal shift
    total_shifts = zeros(size(shifts,2),size(shifts,3));

    for i = 1:size(shifts,3)
        for j = 1:size(shifts,2)
            total_shifts(j,i) = sqrt(shifts(1,j,i)^2+shifts(2,j,i)^2);
        end
    end 

    % Total shifts
    figure
    hold on
    xax = 0:step_size:(step_size*size(total_shifts,2)-step_size);
    for i = 1:size(total_shifts,1)
        plot(xax,total_shifts(i,:),'-','linewidth',2)
    end
    legend({'1', '2', '3', '4', '5', '6', '7'},'Orientation','horizontal','Location','Best');
    box on
    box on
    set(gca,'fontsize',15)
    set(gca,'linewidth',3)
    title(['Lateral Shift Over Axial Distance for' image_params.datafile])
    xlabel('Axial Distance')
    ylabel('Lateral Shift')
    hold off



    % average shift of off axis lenses
    off_ax = reshape(total_shifts(total_shifts>0),[size(total_shifts,1)-1,size(total_shifts,2)]);

    off_ax_av = mean(off_ax,1);
    figure
    plot(xax, off_ax_av,'.','linewidth',2)
    box on
    set(gca,'fontsize',15)
    set(gca,'linewidth',3)
    diff = abs(off_ax_av(end) - off_ax_av(1))/(step_size*size(total_shifts,2));
    title(['Average Off Axis Shift for ' image_params.datafile])%num2str(diff)
    xlabel('Axial Distance')
    ylabel('Lateral Shift')
    hold off

    % Row shifts
    figure
    hold on
    row_shifts = squeeze(shifts(1,:,:));
    for i = 1:size(shifts,2)
        plot(xax,row_shifts(i,:),'.','linewidth',2)
    end
    legend({'1', '2', '3', '4', '5', '6', '7'},'Orientation','horizontal','Location','Best');
    box on
    set(gca,'fontsize',15)
    set(gca,'linewidth',3)
    title(['Row Shift for ' image_params.datafile])
    xlabel('Axial Distance')
    ylabel('Lateral Shift')
    hold off

    % col shifts
    figure
    hold on
    col_shifts = squeeze(shifts(2,:,:));
    for i = 1:size(shifts,2)
        plot(xax,col_shifts(i,:),'.','linewidth',2)
    end
    legend({'1', '2', '3', '4', '5', '6', '7'},'Orientation','horizontal','Location','Best');
    box on
    set(gca,'fontsize',15)
    set(gca,'linewidth',3)
    title(['Column Shift' image_params.datafile])
    xlabel('Axial Distance')
    ylabel('Lateral Shift')
    hold off
end

