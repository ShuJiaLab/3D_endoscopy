function image_params = LDetRect(image_params)
% Lens detection. Image barrel distortion and pixel rectification.

    %% Load image data   
    calpath = image_params.load_paths.cali_path;
    calfile = image_params.load_paths.cali_file;
    cali = double(imread([calpath filesep calfile]));

    %% Diagnostic program to check what might be going wrong   
    % run this if binary image masking isn't working and you want to know
    % why
    calvals = cali(:);
    nonzero = find(calvals);
    calvals = calvals(nonzero);
    if image_params.verbose
        cali_diagnostic
%        image_params.thresh = level;
    end
    
    %% Automatic thresholding or manual thresholding
    if image_params.manual_thresh == 1
        intensity_th = image_params.thresh;
        mask = logical(cali>intensity_th);      
    else
        mask = imbinarize(uint8(cali));
    end
    
    %% Shape detection and size thresholding from binary
    size_th = image_params.size_th;
    shapes = bwconncomp(mask);
    stats = regionprops(shapes, 'Centroid','PixelList',...
        'MajorAxisLength','MinorAxisLength','ConvexImage','boundingbox');
    lengths = cellfun(@length,shapes.PixelIdxList);
    lens_info = stats(lengths>size_th);
    n_lenses = length(lens_info);
    image_params.raw_imsize = size(cali);
    
    %% Generate lens mask image with no padding
    figure
    centers = round(reshape([lens_info.Centroid],[2,7])');
    image_params.circle_radius_nopad = round((mean(mean([lens_info.MajorAxisLength;lens_info.MinorAxisLength]))/2));
    rawim = zeros(size(cali));
    for i=1:n_lenses
        rawim(centers(i,2),centers(i,1)) = 1;   
    end   
    lens_nopad = Generate_RoundBeam(image_params.circle_radius_nopad,2*image_params.circle_radius_nopad+2)./255;
    image_params.lens_mask = conv2(rawim,lens_nopad,'same');
    imagesc(uint8(cali))
    hold on
    viscircles(centers,ones(size(centers,1),1)*image_params.circle_radius_nopad,'Color','r')
    axis off
    colormap(gca,'gray')
    for i=1:n_lenses
        cx = round(lens_info(i).Centroid(1));
        cy = round(lens_info(i).Centroid(2));
        plot(cx,cy,'.r','MarkerSize',20)
    end 
    
    %% Plot lenses with centers
    fig1 = figure;
    imagesc(cali)
    axis off
    colormap(gca,'gray')
    hold on
    for i=1:n_lenses
        cx = round(lens_info(i).Centroid(1));
        cy = round(lens_info(i).Centroid(2));
        plot(cx,cy,'.r','MarkerSize',20)
    end
    

        
    %% Lens Crops - Assign Size
    % the axis of the lens circle based on binary image
    image_params.axes = [lens_info.MajorAxisLength;lens_info.MinorAxisLength];
    % let's pad that somewhat generously and set to an odd pixel length
    padding_bb = 3;
    % We will also want to pad or circle somewhat
    padding_circ = 15;
    % Padding for final BC
    padding_bc = image_params.GRIN_padding;
    
    bb = odd_x(max(max(image_params.axes))+padding_bb);                                          % Find the largest object boundary, add 50 pixels.
    image_params.circle_radius = odd_x(max(max(image_params.axes))/2 + padding_circ);
    % uncorrected lenses will be cropped into this stack
    lenses = zeros(bb,bb,n_lenses);
    % binary images will be placed in this stack
    bin = zeros(bb,bb,n_lenses);
    % the binarized lens images for barrel correction go here
    bin_lenses = zeros(bb,bb,n_lenses);
    % the final bc image will be padding_bc*bb+1 width as determined by the barrel
    % correction code
    bc_width = padding_bc+bb;
    lenses_bc = zeros(bc_width,bc_width,n_lenses);
   
    
    
    %% Make the Crops of Calibration Image
   for i=1:n_lenses
        % bounding box is an odd number, so bb-1/2 gives us the extension 
        % of the box on either side 
        range = (bb-1)/2;
        
        % we round here because the sub pixel information isn't necessarily
        % meaningful or useful
        cx = round(lens_info(i).Centroid(1));
        cy = round(lens_info(i).Centroid(2));
        x = cx - range;
        y = cy - range;
        
        % the extent of our rectangle by upper corner and height/width
        lens_rect(i,:) = [x y bb bb];
        rectangle('Position', lens_rect(i,:),'EdgeColor','r')
        % now we want the actual rows and columns in our box relative to
        % center pixel
        rows = cy-range:cy+range;
        cols = cx-range:cx+range;
        lens_mask_circle = Generate_RoundBeam(image_params.circle_radius,length(rows))./255;
        lenses(:,:,i) = cali(rows,cols).*lens_mask_circle;
        bin(:,:,i) = mask(rows,cols);
        bin_lenses(:,:,i) = cali(rows,cols).*mask(rows,cols);
        lenses_bc(:,:,i) = BarrelDistortionCorrection_V4(lenses(:,:,i),...
            image_params.K1,image_params.K2,image_params.K3,padding_bc);
   end
    
   %% Image mask - no BC
    centers = round(reshape([lens_info.Centroid],[2,7])');
    image_params.circle_radius_nopad = round((mean(mean(image_params.axes))/2));
    rawim = zeros(size(cali));
    for i=1:n_lenses
        rawim(centers(i,2),centers(i,1)) = 1;   
    end   
    lens_nopad = Generate_RoundBeam(image_params.circle_radius_nopad,2*image_params.circle_radius_nopad+2)./255;
    image_params.lens_mask = conv2(rawim,lens_nopad,'same');

    %% Save the barrel correction mask for use with the data
    
    image_params.lens_mask_bc = BarrelDistortionCorrection_V4(lens_mask_circle,image_params.K1,image_params.K2,image_params.K3,padding_bc);

    %% Barrel Distortion Rectification

    % Visualize and re-binarize for lens outline detection
    
    figure
    figpos = get(0,'ScreenSize');
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.2086, 1]);
    ha = tight_subplot(7,3,[.01 .01],[.01 .01],[.01 .01]);
    hold on
    p = 1;
    for i=1:n_lenses
        % show barrel corrected lens
        set(gcf,'CurrentAxes',ha(p))
        imagesc(uint8(lenses_bc(:,:,i)));
        axis off
        xlim([0 1000])
        ylim([0 1000])
        % show gradient image
        [lens_bc_Gmag(:,:,i),lens_bc_Gdir(:,:,i)] = imgradient(lenses_bc(:,:,i),'sobel');
        set(gcf,'CurrentAxes',ha(p+1))
        imagesc(uint8(lens_bc_Gmag(:,:,i)))
        axis off
        xlim([0 1000])
        ylim([0 1000])
        set(gcf,'CurrentAxes',ha(p+2))
        % show binary mask of lens object (this may not always work, as it
        % doesn't for the first step)
        bc_bin(:,:,i) = imbinarize(uint8(lenses_bc(:,:,i)));
        lenses_bc_bin(:,:,i) = lenses_bc(:,:,i).*bc_bin(:,:,i);
        imagesc(bc_bin(:,:,i))
        colormap(gca,'gray')
        axis off
        xlim([0 1000])
        ylim([0 1000])
        p=p+3;
    end
    

  %% bc + raw comparison   

     bc_shapes = struct();
     bc_lens_bin = zeros(size(lenses_bc_bin));
     for i=1:n_lenses

        % generate bc_lens binary images
        bc_shapes(i).shape = bwconncomp(bc_bin(:,:,i));
        bc_obj = bc_shapes(i).shape.PixelIdxList;
        bc_obj_lens = cellfun(@length,bc_obj);
        bc_lens = bc_obj{bc_obj_lens == max(bc_obj_lens)};
        bc_lens_idx = zeros(size(bc_bin,1),size(bc_bin,2));
        bc_lens_idx(bc_lens) = 1;
        bc_lens_bin(:,:,i) = bc_lens_idx;
        bc_shapes(i).stats = regionprops(bc_lens_bin(:,:,i), 'Centroid','PixelList',...
    'MajorAxisLength','MinorAxisLength','ConvexImage','boundingbox');
        bc_shapes(i).centroid = bc_shapes(i).stats.Centroid;
        bc_shapes(i).bb = bc_shapes(i).stats.BoundingBox;
        if image_params.verbose == 1
            p = 1;
            figure
            ha = tight_subplot(2,2,[.001 .001],[.01 .01],[.01 .01]);
            % plot raw intensity image
            set(gcf,'CurrentAxes',ha(p))
            axis off
            imshow(uint8(lenses(:,:,i)));
            colormap(gca,'hot')
            p = p+1;
            % plot raw binary
            set(gcf,'CurrentAxes',ha(p))
            axis off
            imshow(bin(:,:,i));
            p = p+1;        
            % plot bc corrected intensity
            set(gcf,'CurrentAxes',ha(p))
            axis off
            imshow(uint8(lenses_bc(:,:,i)))
            colormap(gca,'hot')
            p = p+1;
            % plot bc corrected binary
            set(gcf,'CurrentAxes',ha(p))
            axis off
            imshow(bc_lens_bin(:,:,i))
            p = p+1;     
         end 
     end
        
    %% Centroid rectification
%     
     raw_centroids = reshape([lens_info.Centroid],2,7)';
%     imsize = size(cali);    
     bc_size = size(bc_bin(:,:,1));
     bc_r = (bc_size(1)-1)/2;
%     r_c = round(raw_centroids(:,2));
%     c_c = round(raw_centroids(:,1));
%     row_max = imsize(1);
%     col_max = imsize(2);
%     
%     % left col clearance
%     rcc = abs(bc_r-min(min(c_c)));
%     % right col clearance
%     lcc = abs(bc_r-(col_max - max(max(c_c))));
%     % top row clearance
%     trc = abs(bc_r-min(min(r_c)));
%     % bottom row clearance
%     brc = abs(bc_r - (row_max - max(max(r_c)))); 
%     
%     bc_padding = max([rcc,lcc,trc,brc]);
%     bc_padding = ceil(bc_padding)+5;
    
    largest_lens_spacing = max(max(abs(raw_centroids(4,:)-raw_centroids)));
    fullimsize = odd_x(2*(largest_lens_spacing+bc_r));
    image_params.bc_fullimsize = [fullimsize, fullimsize];
    
    new_center = [ceil(fullimsize/2),ceil(fullimsize/2)];
    bc_fullsize_centroids = abs(new_center - ((raw_centroids(4,:)-raw_centroids)));
    
    % visualize bc and raw image overlaid for sanity check
    bc_cal = zeros(image_params.bc_fullimsize(1),image_params.bc_fullimsize(2),3);
    for i =1:n_lenses
        c = odd_x(bc_fullsize_centroids(i,:));
        image_params.bc_fullsize_centroids(i,:) = c;
        raw_rows = c(2) - range : c(2) + range;
        raw_cols = c(1) - range : c(1) + range;
        bc_row = c(2) - bc_r : c(2) + bc_r;
        bc_cols = c(1) - bc_r : c(1) + bc_r;

        bc_cal(c(2)-20:c(2)+20,c(1)-20:c(1)+20,1) = 255;
        bc_cal(raw_rows,raw_cols,2) = lenses(:,:,i);
        bc_cal(bc_row,bc_cols,3) = lenses_bc(:,:,i);
    end    

    figure
    imshow(uint8(cat(3,bc_cal)));
    
    
    %% Populate Output 
    
    for i = 1:n_lenses
        lens_info(i).lenses_bc = lenses_bc(:,:,i);
        lens_info(i).lenses_bc_bin_only = bc_bin(:,:,i);
        lens_info(i).lenses_bc_bin_int = lenses_bc_bin(:,:,i);     
    end
    
    image_params.lens_info = lens_info;
    image_params.raw_lens_crop = bb;
end