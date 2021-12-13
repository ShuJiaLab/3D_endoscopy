%% A snippet from RayOpticsRecon
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