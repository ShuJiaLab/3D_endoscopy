function image_params = Shift_Caclulation(image_params,psf_overlay) 
% Create three new vector quantities
% shift_dir_vec - direction along which to shift
% starting_pos - position for each image to start at
% step_shift - shift for each PSF step

%% Create a pictoral representation of the PSF and of the shift centers to get a sense of what's going on

    % A circle for visualization
    GRIN_BC_R = ceil(size(image_params.lens_info(1).lenses_bc_bin_only,1)/2)
    Output = Generate_RoundBeam( GRIN_BC_R , GRIN_BC_R*2+1 )/255; 
    Output_2 = Generate_RoundBeam( GRIN_BC_R-5 , GRIN_BC_R*2+1 )/255; 
    Mask_Circle          = Output - Output_2;
    
    % PSF Visualization
    %image_mask = zeros
    image_mask = zeros(image_params.bc_fullimsize);
    %centers = reshape(round([image_params.lens_info.Centroid]),[2,7])';
    if image_params.manual_lens_fit == 1
        %%%%% MANUAL SHIFT DATA %%%%%
        % Note
        %centers(:,1) = center_row;
        %centers(:,2) = center_col;  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        manualshift
        image_params.manual_centers = centers;
        image_params.bc_manual_centers = bc_centers;

        centers = image_params.bc_manual_centers;
    else
        centers = image_params.bc_fullsize_centroids;
    end
    
    for i=1:length(centers)
        image_mask(centers(i,2),centers(i,1)) = 1;
    end
    lenses = 2*conv2(image_mask , Mask_Circle , 'same')*255;

    
    lens_psf = zeros(size(image_mask,1),size(image_mask,2),3);
    lens_psf(:,:,1) = lenses;
    lens_psf(:,:,2) = psf_overlay;
    lens_psf(:,:,3) = 2*conv2(image_mask, Generate_RoundBeam(5,11),'same')*255;
    
    imshow(uint8(lens_psf))
    
    psf_length = size(image_params.Positionsxy,3);
    lens_num = length(image_params.lens_info);
    
    
    
    shifts = zeros(size(image_params.Positionsxy));
    shifts_normalized = zeros(size(image_params.Positionsxy));
    for i = 1:psf_length % we will start with the second step of the psf and compute relative shift for each point
        
        for j = 1:lens_num 
            shifts(1,j,i) = image_params.Positionsxy(1,j,i) - image_params.bc_fullsize_centroids(j,1); % columns
            shifts(2,j,i) = image_params.Positionsxy(2,j,i) - image_params.bc_fullsize_centroids(j,2); % rows
        end
        cs = shifts(:,:,i);
        shifts_normalized(:,:,i) = cs - cs(:,4);
    end    
    
    image_params.shifts = shifts;
    image_params.shifts_normalized = shifts_normalized;
end
