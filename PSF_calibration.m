function [image_params,psf_overlay] = PSF_calibration(image_params)
    
    %%% PSF INPUTS %%%
    psf_folder = dir(image_params.load_paths.psf_path);
    %psf_folder = psf_folder([psf_folder.bytes]>0);
    psf_folder = psf_folder(~[psf_folder(:).isdir]);
    
    
    %%% LOAD OPTIONAL CALIBRATION DATA %%%
    if image_params.PSF_cal_in
        load([image_params.load_paths.psf_cal_in_file])
        cal_positions = Positionsxy;
    end
    
    %%% PSF PARAMETERS %%% 
    psf_depth = length(psf_folder);
    image_params.psf_depth = psf_depth;
    num_GRIN = length(image_params.lens_info);
    di_R = (image_params.raw_lens_crop - 1)/2;
    Pixel_D = image_params.Pixel_D;
    
    %%% GAUSSIAN FITTING %%%
    x                  = -24:1:24;
    x                  = x*Pixel_D;
    y                  = x;
    [Xg Yg] = meshgrid(x ,y);
    Pixel_Num                = 1.5;%sqrt( ( (Index - 60)^2 + 1280 )/400 );   % --- > Does this gaussian width still make sense?
    Q                        = Pixel_Num*Pixel_D;
    Gaussian                 = 1/(Q*sqrt(2*pi)) * exp( -1/2 * ( ( Xg/Q ).^2 + ( Yg/Q ).^2 ) );
    Gaussian                 = Gaussian/max(Gaussian(:))*255;

    GRIN = Generate_RoundBeam( di_R, (di_R)*2+1 )/255; % GRIN outline
    Positionsxy = zeros(2,num_GRIN);
    
    figure
    
    % PSF Overlay will be the size of the barrel corrected image
    psf_overlay = zeros(image_params.bc_fullimsize);
    bc_R = (size(image_params.lens_info(1).lenses_bc,1)-1)/2;
    for i=1:psf_depth % For each layer of the PSF
        psf_layer = im2double(imread([psf_folder(i).folder filesep psf_folder(i).name]));
        [l3,~,~,t3] = islink([psf_folder(i).folder filesep psf_folder(i).name]);
        if l3
            image_params.load_paths.psf_file_link{i} = t3;
        end
        gauss_crop=zeros(2*image_params.PSF_gaussbox+1,2*image_params.PSF_gaussbox+1,num_GRIN);  
        positions_output = zeros(image_params.bc_fullimsize);
        for a=1:num_GRIN % For each GRIN lens
            GRINcenter = zeros(image_params.raw_imsize);
            c_c = round(image_params.lens_info(a).Centroid(1));
            c_r = round(image_params.lens_info(a).Centroid(2));
            c_c_bc = round(image_params.bc_fullsize_centroids(a,1));
            c_r_bc = round(image_params.bc_fullsize_centroids(a,2));
            GRINcenter(c_r,c_c) = 1;
            GRINbounds = conv2(GRINcenter, GRIN, 'same')*1;
            %imshow(GRINbounds)
            GRIN_psf = psf_layer.*GRINbounds;
            GRIN_psf_bc = zeros(image_params.bc_fullimsize);
           
            %%% BARREL DISTORTION CORRECTION %%%
            GRIN_psf_bc(c_r_bc-bc_R:c_r_bc+bc_R,c_c_bc-bc_R:c_c_bc+bc_R)...
                = BarrelDistortionCorrection_V4( GRIN_psf(c_r-di_R:c_r+di_R,c_c-di_R:c_c+di_R),...
                image_params.K1,image_params.K2,image_params.K3,image_params.GRIN_padding);
            
            disp(i)
 
            %%% PSF THRESHOLDING %%%
            PSFvalues = GRIN_psf_bc(:);
            nonzero = find(PSFvalues);
            PSFvalues = PSFvalues(nonzero);
            if (size(GRIN_psf_bc,1) ~= size(positions_output,1)|size(GRIN_psf_bc,2) ~= size(positions_output,2))
                [r, c] = size(positions_output);
                GRIN_psf_bc = imresize(GRIN_psf_bc, [r, c]);
            end
            % % Removing thresholding temporarily
%             %med = median(PSFvalues);
%             %thresh = 2.75*med; % may need to change to higher if high intensity stray light, 3*med didn't work for all slides in one run
%             %PSF_thresh = GRIN_psf_bc > thresh; %creates an array with 1 where the intensity is above the threshold
%             %mask = zeros(size(GRIN_psf_bc));
%             %mask(PSF_thresh) = imbinarize(GRIN_psf_bc(PSF_thresh));
            % % End thresholding
            
            mask = imbinarize(GRIN_psf_bc);
            shapes = bwconncomp(mask);
            stats = regionprops(shapes,'Centroid');
            objects = shapes.PixelIdxList;
            lengths = cellfun(@length,objects);
            [M, I] = max(lengths);
            temp2 = {stats.Centroid};
            CenterPixel = temp2(I);
            
            % Center pixel of object
            if ~image_params.PSF_cal_in
                CenterPixelArray = cell2mat(CenterPixel);
            elseif image_params.PSF_cal_in
            % Replace center pixel with calibration pixels for low SNR PSFs
                CenterPixelArray = [cal_positions(1,a,i) cal_positions(2,a,i)];
            end
            gauss_crop(:,:,a) = GRIN_psf_bc(CenterPixelArray(1,2)-image_params.PSF_gaussbox:CenterPixelArray(1,2)+image_params.PSF_gaussbox , CenterPixelArray(1,1)-image_params.PSF_gaussbox:CenterPixelArray(1,1)+image_params.PSF_gaussbox);

            x0 = [1, 0, 50, 0, 50, 0]; %initial guess params
            [~, fitpara] = GaussianFitting_2D(gauss_crop(:,:,a), x0, 100);
            PSF_r = round(CenterPixelArray(1,2))+round(fitpara(1,4));
            PSF_c = round(CenterPixelArray(1,1))+round(fitpara(1,2));
            
            PSF_temp = zeros(size(GRIN_psf_bc));
            
            %%% RECORD OUTPUT POSITIONS %%%
            Positionsxy(1, a, i) = PSF_c;
            Positionsxy(2, a, i) = PSF_r;

            image_params.fitpara(a,i) = {fitpara};
            
            % Record output positions pictorally
            PSF_temp(PSF_r, PSF_c) = 1;    
            
            PSF_temp = conv2(PSF_temp,Gaussian,'same');
            
            if size(PSF_temp) ~= size(positions_output)
                [r, c] = size(positions_output);
                PSF_temp = PSF_temp(1:r,1:c);
            end
           
            positions_output = positions_output + PSF_temp;
            %Gauss Center Overlay Image making
            PSF_center = zeros(size(GRIN_psf_bc));
            PSFCenterCol=round(CenterPixelArray(1));
            PSFCenterRow=round(CenterPixelArray(2));
            PSF_center(PSFCenterRow-1:PSFCenterRow+1,PSFCenterCol-1:PSFCenterCol+1) = 255;
            PSFwGaussCenter(:,:,1) = PSF_center;
            PSFwGaussCenter(:,:,2)      = GRIN_psf_bc;
            PSFwGaussCenter(:,:,3)      = 0;
            imwrite(PSFwGaussCenter(PSFCenterRow-30:PSFCenterRow+30,PSFCenterCol-30:PSFCenterCol+30,:),[image_params.save_paths.psf_save filesep 'GaussCenterOverlay.tif'],'compression','none','writemode','append');
        end


        
        %imwrite(uint8(positions_output),[image_params.save_paths.psf_save filesep 'positions_matrix.tif'],'compression','none','writemode','append');
        psf_overlay = psf_overlay + positions_output;
        
%         if image_params.long_output == 1         
%             disp_psf_layer(gauss_crop,psf_save)
%         end      
    end
    
    image_params.Positionsxy = Positionsxy;
    
    save([image_params.save_paths.psf_save filesep 'positionsxy.mat'],'Positionsxy')
    imshow(psf_overlay);
    imwrite(uint8(psf_overlay),[image_params.save_paths.psf_save filesep 'psf_overlay.tif'],'compression','none','writemode','append');
   
end