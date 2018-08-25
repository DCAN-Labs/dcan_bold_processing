function format_string = format_generator(skip_keep_vector)
% This function assumes 1 = remove or skip and 0 = keep. It also assumes
% the pre-skip frame and 2 post-skip frames are already in the
% skip_keep_vector.
if isstr(skip_keep_vector)
    disp('You are silly, because you put a string, when I want a vector')
    return
end

if size(skip_keep_vector,1)>1
    skip_keep_vector = skip_keep_vector'; 
end 

indices = [0 find(diff(skip_keep_vector)~=0)];

format_string = [];

for i = 1:length(indices)
    if i ~= length(indices)
        temp = skip_keep_vector(indices(i)+1:indices(i+1));
        format_string = cat(2,format_string,num2str(length(temp)));
    else
        temp = skip_keep_vector(indices(i)+1:length(skip_keep_vector));
        format_string = cat(2,format_string,num2str(length(temp)));
    end
    if sum(temp)~=0
        format_string = cat(2,format_string,'x');
    else
        format_string = cat(2,format_string,'+');
    end
end
