%% Written by Eric Earl, 10/15/2014
function subject_motion_numbers_TXT_parse_v2_octave(movement_files,output_dir)
% Motion numbers parsing Nipype wrapper pulling in list of motion files
% Outputs motion_numbers to "motion_numbers.mat" file
% movement_files = {'motion_numbers0.txt', 'motion_numbers1.txt', 'motion_numbers2.txt'};
% movement_files must be "motion_numbers.txt" files from CYAPP pipelines in a MATLAB cell of comma-separated strings

run_count = size(movement_files,2);

% Initialize the empty vectors for storing the parsed motion_numbers
motion_numbers.diff_x = [];
motion_numbers.diff_y = [];
motion_numbers.diff_z = [];
motion_numbers.diff_arc_len_x = [];
motion_numbers.diff_arc_len_y = [];
motion_numbers.diff_arc_len_z = [];
motion_numbers.FD = [];
motion_numbers.DVAR_pre_reg = [];
motion_numbers.DVAR_post_reg = [];
motion_numbers.DVAR_post_all = [];

%%
for i = 1:run_count
    movement_file  = movement_files{i};
    
    % Read movement files for backwards difference "motion" numbers
    % Pad vectors with a false frame (difference of frame "0" and 1)
    % before the first numbers (the difference of frame 1 and 2)
    switch movement_file(end-3:end)
        case '.txt'
            fid = fopen(movement_file);
            motion_cell = textscan(fid, '%*s %*s %*s %f %f %f %f %f %f %f %f %f %f','HeaderLines',1);% Read FD values from motion numbers files
            fclose(fid);
            diff_x = cat(1,0,motion_cell{1});
            diff_y = cat(1,0,motion_cell{2});
            diff_z = cat(1,0,motion_cell{3});
            diff_arc_len_x = cat(1,0,motion_cell{4});
            diff_arc_len_y = cat(1,0,motion_cell{5});
            diff_arc_len_z = cat(1,0,motion_cell{6});
            FD = cat(1,0,motion_cell{7});
            DVAR_pre_reg = cat(1,0,motion_cell{8});
            DVAR_post_reg = cat(1,0,motion_cell{9});
            DVAR_post_all = cat(1,0,motion_cell{10});
        otherwise
            disp('Input is not a .txt file.  Exiting.')
            return
    end
    
    motion_numbers.diff_x = cat(1,motion_numbers.diff_x,diff_x);
    motion_numbers.diff_y = cat(1,motion_numbers.diff_y,diff_y);
    motion_numbers.diff_z = cat(1,motion_numbers.diff_z,diff_z);
    motion_numbers.diff_arc_len_x = cat(1,motion_numbers.diff_arc_len_x,diff_arc_len_x);
    motion_numbers.diff_arc_len_y = cat(1,motion_numbers.diff_arc_len_y,diff_arc_len_y);
    motion_numbers.diff_arc_len_z = cat(1,motion_numbers.diff_arc_len_z,diff_arc_len_z);
    motion_numbers.FD = cat(1,motion_numbers.FD,FD);
    motion_numbers.DVAR_pre_reg = cat(1,motion_numbers.DVAR_pre_reg,DVAR_pre_reg);
    motion_numbers.DVAR_post_reg = cat(1,motion_numbers.DVAR_post_reg,DVAR_post_reg);
    motion_numbers.DVAR_post_all = cat(1,motion_numbers.DVAR_post_all,DVAR_post_all);
    
end

%%
save([output_dir filesep 'motion_numbers.mat'],'motion_numbers','-v7');
