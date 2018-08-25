function [scanner_ss,count_scanner,count_processed] = sweetspot(scanner,processed,th_index)

%load(name)

%Define range,number of FD thresholds
FD = [0.01:0.01:0.5];
nth = length(FD);
th = FD(th_index); %threshold for FD_processed to determine remaining time
n_vol = length(scanner);%Determine paramters of input varibles

%preallocate variables
count_scanner = zeros(1,nth);
count_processed = zeros(1,nth);

for i=1:nth; 
    count_scanner(i) = sum(scanner<FD(i)); 
end
    
for i=1:nth; 
    count_processed(i) = sum(processed<FD(i)); 
end

%get the remaining frame count for the processed
remaining_processed = count_processed(th_index);

%see at what threshold the scanner fd has the same/similar frame count
index = min(find(count_scanner>remaining_processed));
scanner_ss = FD(index);

%find y value (cumulative time) to plot a line at the sweetspot point
%y_val = count_scanner(index)

%plot figure (cumulative distribution...sort of)
%figure;
%plot(1:50,count_scanner,1:50,count_processed)
%plot(1:nth,count_scanner,1:nth,count_processed)
%grid on 
%grid minor
%title(name)
%legend('Scanner Remaning', 'Processed Remaining', 'location', 'SE')
%line([0 nth], [y_val y_val])


