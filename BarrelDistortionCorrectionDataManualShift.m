
function image_params = BarrelDistortionCorrectionDataManualShift(image_params,data)
   
    %parameters
    num_GRIN = length(image_params.lens_info)

    data_bc = zeros(size(image_params.bc_fullimsize));
    bb = image_params.raw_lens_crop;
    range = (bb-1)/2;
    
    f1 = figure;
    imshow(uint8(data))
    
    f2 = figure;
    figpos = get(0,'ScreenSize');
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.2086, 1]);
    ha = tight_subplot(7,3,[.01 .01],[.01 .01],[.01 .01]);
    p = 1;    
    for a=1:num_GRIN
        cx = image_params.manual_centers(a,1);
        cy = image_params.manual_centers(a,2);
        x = cx - range;
        y = cy - range;
        set(0,'currentfigure',f1)
        hold on
        plot(cx,cy,'.r','MarkerSize',20)

        % the extent of our rectangle by upper corner and height/width
        lens_rect(a,:) = [x y bb bb];
        rectangle('Position', lens_rect(a,:),'EdgeColor','r')
        % now we want the actual rows and columns in our box relative to
        % center pixel
        rows = cy-range:cy+range;
        cols = cx-range:cx+range; 
        %% Take lens lens crop
        lens_i(:,:,a) = data(rows,cols);
        
        % Visualize
        set(0,'currentfigure',f2)
        hold on
        set(gcf,'CurrentAxes',ha(p))
        imagesc(uint8(lens_i(:,:,a)));
        axis off
        xlim([0 1000])
        ylim([0 1000])
        
        %% Barrel Correct
        bc_lens_i(:,:,a)...
        = BarrelDistortionCorrection_V4(lens_i(:,:,a),....
        image_params.K1,image_params.K2,image_params.K3,...
        image_params.GRIN_padding);
    
        % Visualize
        set(gcf,'CurrentAxes',ha(p+1))
        imagesc(uint8(bc_lens_i(:,:,a)));
        axis off
        xlim([0 1000])
        ylim([0 1000])    

        %% Binarize
        bc_lens_i_thresh(:,:,a) = ...
            bc_lens_i(:,:,a) .* image_params.lens_mask_bc;
        
        % Visualize
        set(gcf,'CurrentAxes',ha(p+2))
        imagesc(uint8(bc_lens_i_thresh(:,:,a)));
        axis off
        xlim([0 1000])
        ylim([0 1000])         
         
        p = p+3;
    end 

    image_params.bc_data_nothresh = bc_lens_i;
    image_params.bc_data_circthresh = bc_lens_i_thresh;

end