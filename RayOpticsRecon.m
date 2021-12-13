 function [Recon, image_params] = RayOpticsRecon(image_params)
% For first layer of reconstruction, add all centers together
% For second layer, displace each center by delta d, to the opposite side
% of the reconstruction center as the lens is in physical space. Delta d is
% given by displacements of the system response function
    
    %% Initialize Reconstruction Parameters
    if image_params.normshifts == 1
        shifts = image_params.shifts_normalized;
    elseif image_params.normshifts == 0
        shifts = image_params.shifts;
    end
    
    psf_length = size(image_params.Positionsxy,3);  
    num_lenses = size(image_params.Positionsxy,2); 
   
        
    % Diameter of one subimage
    subim_D = size(image_params.bc_data_nothresh,1)-1;
    
    image_params.max_row_shift = max(max(max(abs(shifts(1,:,:)))));
    image_params.max_col_shift = max(max(max(abs(shifts(2,:,:)))));
    subim_R = round(subim_D/2);
    %padding = image_params.GRIN_padding
    m_s_r = image_params.max_row_shift;
    m_s_c = image_params.max_col_shift;
    m_s = max(m_s_r, m_s_c);
    psf_length = size(image_params.Positionsxy,3);  
    num_lenses = size(image_params.Positionsxy,2);
    Recon_rows = odd_x(2*(m_s+subim_R));
    Recon_cols = odd_x(2*(m_s+subim_R));
    num_lenses = size(image_params.Positionsxy,2); 

    Recon = zeros(Recon_rows, Recon_cols, psf_length);

    % Sizing


    %subim_R = round(subim_D/2);
    
    % Image center - assumes 7 lenses
    % center of reconstruction
    cc = ceil(Recon_cols/2);
    % reconstruction edge
    re = cc - subim_R;
    
    %rtr = cc - GRIN_R - padding; % recon top row
    %rec = cc - GRIN_R - padding; % recon edge column
    
    if image_params.MAG_cal == 1
        
        load(image_params.load_paths.mag_pos_file)
        M = rgb_shifts.avg_m_allcorlor_f;
        mags = M(image_params.recon_range + M.o);
        shifts = mags*image_params.lens_pitch./(image_params.Pixel_D);
        
        tic
        for i = 1:psf_length
            for j=1:num_lenses
               lcr = image_params.bc_fullsize_centroids(j,1); % lens center row
               lcc = image_params.bc_fullsize_centroids(j,2); % lens center column
               Recon(re - shifts(2,j,i):re - shifts(2,j,i)+subim_D, re - shifts(1,j,i):re - shifts(1,j,i)+subim_D,i) = ...
                Recon(re - shifts(2,j,i):re - shifts(2,j,i)+subim_D, re - shifts(1,j,i):re - shifts(1,j,i)+subim_D,i)+...        
                image_params.bc_data_nothresh(:,:,j);
            end
            Recon(:,:,i) = (Recon(:,:,i)./num_lenses);
        end    
        toc
        
    else
    %% Shift and add them to the reconstruction
        tic
        for i = 1:psf_length
            for j=1:num_lenses
               lcr = image_params.bc_fullsize_centroids(j,1); % lens center row
               lcc = image_params.bc_fullsize_centroids(j,2); % lens center column
               Recon(re - shifts(2,j,i):re - shifts(2,j,i)+subim_D, re - shifts(1,j,i):re - shifts(1,j,i)+subim_D,i) = ...
                Recon(re - shifts(2,j,i):re - shifts(2,j,i)+subim_D, re - shifts(1,j,i):re - shifts(1,j,i)+subim_D,i)+...        
                image_params.bc_data_nothresh(:,:,j);
            end
            Recon(:,:,i) = (Recon(:,:,i)./num_lenses);
        end    
        toc
    end
%     h = imshow(uint8(Recon(:,:,1))); set(h,'erasemode','xor');
%     for i=2:psf_length
%         set(h,'cdata',Recon(:,:,i))
%         pause
%     end
    
    t = datetime;
    t.Format = 'MM_dd_uuuu_hh_mm_ss';
    
    for i = 1:psf_length
        imwrite(uint8(Recon(:,:,i)),...
            strcat(image_params.save_paths.recon_save, filesep, 'recon_', image_params.datafile, '_', string(t), '.tif'),'compression','none','writemode','append');    
    end
    
    save(strcat(image_params.save_paths.recon_save, filesep, 'image_params_', image_params.datafile, '_', string(t), '.mat'),'image_params')
    
 end
 
 %image_params.max_col_shift = max(max(max(abs(shifts(2,:,:)))));

     %image_params.max_row_shift = max(max(max(abs(shifts(1,:,:)))));
         %GRIN_R = image_params.GRIN_R_S;
    %padding = image_params.GRIN_padding
    %m_s_r = image_params.max_row_shift;
    %m_s_c = image_params.max_col_shift;
    %m_s = max(m_s_r, m_s_c);
    %     rcr = round(size(Recon,1)/2); % recon center row
%     rcc = round(size(Recon,2)/2); % recon center column
%     
%     rtr = rcr - GRIN_R - padding; % recon top row
%     rec = rcc - GRIN_R - padding; % recon edge column
%     
%     subim_D = 2*GRIN_R + 2*padding;
%     subim_R = subim_D/2;
    

    %intersect_map = Intersect_Calculation(GRIN_R, rcr, rcc, shifts, Recon_rows, Recon_cols);