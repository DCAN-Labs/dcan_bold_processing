function sorted_path=sorter(defined_path,taskname)

path_contents=dir(defined_path);
nn=size(path_contents,1);
ix_ls=zeros(nn,1);
sorted_path=cell(nn,1);

for ii=1:nn
    filename=[path_contents(ii).folder filesep path_contents(ii).name];
    [tasknumstarts,tasknumends] = regexp(filename,[taskname '[0-9]{2}']);
    ix_ls(ii)=str2num(filename((tasknumends(end)-1):tasknumends(end)));
end

[cc, ix]=sort(ix_ls);

for ii=1:nn
    path_contents(ix(ii)).folder;
    sorted_path{ii}=[path_contents(ix(ii)).folder filesep path_contents(ix(ii)).name];
end
