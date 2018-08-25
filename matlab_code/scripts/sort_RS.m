function sorted_path=sort_RS(defined_path)

path_contents=dir(defined_path);
nn=size(path_contents,1);
ix_ls=zeros(nn,1);
sorted_path=cell(nn,1);

for ii=1:nn
    filename=[path_contents(ii).folder filesep path_contents(ii).name];
    [restnumstarts,restnumends] = regexp(filename,'REST[0-9]+');
    ix_ls(ii)=str2num(filename((restnumstarts(end)+4):restnumends(end)));
end

[cc, ix]=sort(ix_ls);

for ii=1:nn
    path_contents(ix(ii)).folder;
    sorted_path{ii}=[path_contents(ix(ii)).folder filesep path_contents(ix(ii)).name];
end

