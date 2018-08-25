function raw_motion = extract_MoCo(dicomdir)
%% This script parses the movement information from MoCo series
% NOTICE: The raw motion output will be in the native ordering format
% (y,x,z,y_rot,x_rot,z_rot)


%dicomdir = '/group_shares/FAIR_ADHD/CYA/sorted/ADHD-HumanYouth-OHSU/11778-1/20150211-SIEMENS_TrioTim-Nagel_K_Study/1755-6-REST3-MoCoSeries';
N = length(dir(dicomdir))-2;
filenames = 0;
current_files = {};
raw_motion = zeros(N,6);

%% While loop iterates over the total number of dicoms
while length(filenames) < N
    % list directory, get new files
    dicom_files = dir(dicomdir);
    filenames = {dicom_files(3:end).name};
    new_files = setdiff(filenames,current_files);
    if ~isempty(new_files)
        % deal with the new files
        Nfiles = size(new_files,2);
        disp(['Reading ' num2str(Nfiles) ' new DICOM files...']);
        
        for i=1:Nfiles
            dcmfile = [dicomdir filesep new_files{i}];
            [~,~,ext] = fileparts(dcmfile);
           
            % check to see if the new files are DICOM (or ima)
            if ( strcmp(ext,'.dcm') || strcmp(ext,'.IMA') || strcmp(ext,'.ima') )
                dcmfile = [dicomdir filesep new_files{i}];
               
                % use image processing toolbox function "dicominfo" to get header info
                info = dicominfo(dcmfile);
               
                % if there is Siemens MoCo info in the header, there will be an "ImageComments" field
                if ( isfield(info,'ImageComments') )
                    frame = info.AcquisitionNumber;
                    if ( strncmp(info.ImageComments,'Motion:',7) )
                        raw_motion(frame,:) = sscanf(info.ImageComments,'Motion: %f,%f,%f,%f,%f,%f');
                    end
                end
            end
        end
        current_files = filenames;
       
    else
        disp('No new DICOMs.');
        
    end
   
end
 