
function image_params = BarrelDistortionCorrectionData(image_params,data)
   
    %parameters
    num_GRIN = length(image_params.lens_info)

    data_bc = zeros(size(image_params.bc_fullimsize));
    bb = image_params.raw_lens_crop;
    range = (bb-1)/2;
    
    % find lens position from centroids
    c_raw_im = reshape([image_params.lens_info(:).Centroid],[2,7])'
    image_params.positions = find_positions_from_centroids(c_raw_im);
    
    f1 = figure;
    imshow(uint8(data))
    
    f2 = figure;
    figpos = get(0,'ScreenSize');
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.2086, 1]);
    ha = tight_subplot(7,3,[.01 .01],[.01 .01],[.01 .01]);
    p = 1;    
    for a=1:num_GRIN
        cx = round(image_params.lens_info(a).Centroid(1));
        cy = round(image_params.lens_info(a).Centroid(2));
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

    % Make a barrel distortion corrected data image
    
%     bc_sidelength = size(image_params.bc_data_circthresh,1)
%     bc_width = (bc_sidelength-1)/2;
%     
%     bc_fullim = zeros(size(image_params.bc_fullimsize));
%     
%     pos = image_params.positions.lens_positions;
%     lens_position = fieldnames(pos)
%     %centroids_bc = image_params.bc_fullsize_centroids;
%     for i=1:size(image_params.bc_data_circthresh,3)
%         centroid = pos.(lens_position{i})
%         match_idx(i) = find(ismember(c_raw_im,centroid),1)
%         bc_row = image_params.bc_fullsize_centroids(match_idx(i),2);
%         bc_col = image_params.bc_fullsize_centroids(match_idx(i),1);
%         image_params.positions.lens_positions.(lens_position{i}) = [bc_col, bc_row]; 
%         corr_image = image_params.bc_data_circthresh(:,:,match_idx(i));
%         rows = bc_row-bc_width:bc_row+bc_width;
%         cols = bc_col-bc_width:bc_col+bc_width;
%         bc_fullim(rows,cols) = corr_image;
%     end
%     
%     figure
%     
%     imshow(uint8(bc_fullim))
    
    image_params.bc_data_nothresh = bc_lens_i;
    image_params.bc_data_circthresh = bc_lens_i_thresh;
    %image_params.bc_fullim = bc_fullim;
    %image_params.positions.lens_positions.match_idx = match_idx;

end