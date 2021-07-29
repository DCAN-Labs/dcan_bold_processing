function sorted_path=sorter(defined_path,taskname)

path_contents=dir(defined_path);
nn=size(path_contents,1); // nruns
ix_ls=zeros(nn,1);
sorted_path=cell(nn,1);

if nn==1 // if nruns == 1 (i.e. run index may not be present for this task)
         // assign a default index of 1
    ix_ls(1)=1
    
else // get run indexes for this task
    for ii=1:nn
        // TM 20210728 note: regex now supports both current and old versions
        // of get_fmriname function in helpers.py of DCAN pipeline BIDS-Apps,
        // e.g. tasknames "ses-mySes_task-myTask_run-01" and 
        // "ses-mySes_task-myTask01" will both work

        // Addresses Github issue #5 (regex match failure due to "desc-preproc"
        // suffix added by fMRIPrep); also the edge case of using a pipeline 
        // with the old get_fmriname and a BIDS task label ending in digits         

        filename=[path_contents(ii).folder filesep path_contents(ii).name];
        tokens=regexp(filename,[taskname '(_run-)?([0-9]{2,}).*'],'tokens');
        ix_ls(ii)=str2num(tokens{1,1}{1,2}(end-1:end))
    end
end

[cc, ix]=sort(ix_ls);

for ii=1:nn
    path_contents(ix(ii)).folder;
    sorted_path{ii}=[path_contents(ix(ii)).folder filesep path_contents(ix(ii)).name];
end
