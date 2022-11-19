function label_map = get_label_map(moby_anatomy_prefix, moby_lesion_prefix, N, frame_i)

    % Load MOBY phantom and spherical lesion phantom
    fname_body = fullfile(moby_anatomy_prefix, ...
        ['moby_act_', num2str(frame_i), '.bin']);
    fname_lesn = fullfile(moby_lesion_prefix, ...
        ['lesion_act_', num2str(frame_i), '.bin']);
    
    fid = fopen(fname_body, 'rb'); phan = fread(fid, 'float'); fclose(fid);
    fid = fopen(fname_lesn, 'rb'); lesn = fread(fid, 'float'); fclose(fid);
    
    % Insert the spherical lesion into the phantom
    phan(lesn ~= 0) = lesn(lesn ~= 0);
    label_map = reshape(phan, N);
end

