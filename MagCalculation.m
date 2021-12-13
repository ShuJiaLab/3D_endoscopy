function image_params = MagCalculation(image_params)

    mag = struct()
    fit_params = struct()
    % lens size
    lens_D = image_params.circle_radius*2;
    obj_pix_size = 1/lens_D;
    im_pix_size = image_params.Pixel_D;
    
    PlotShiftsOverDistance

    zsteps = length(off_ax_av);
    zstepsize = image_params.step_size;
    zrange = zstepsize:zstepsize:zsteps*zstepsize;
    pitch = image_params.lens_pitch;
    
    magnification = (off_ax_av.*im_pix_size)/image_params.lens_pitch;

    plot(zrange,magnification)
    
    % Fit
    x = zrange';
    y = magnification';
    ft = fittype( 'magfit(x,a,o)');
    f = fit(x,y,ft,'StartPoint',[15,3]);
    
    fit_params.x = x;
    fit_params.y = y;
    fit_params.f = f;
    
    mag.off_ax_shifts = off_ax;
    mag.off_ax_av_shifts = off_ax_av;
    
    image_params.mag = mag;
    image_params.fit = fit_params;
end