New_fullV(1:length(New_S3{1}(:,37))) = New_S3{1}(:,37);
Old_fullV(1:length(Old_S3{1}(:,37))) = Old_S3{1}(:,37);


for i = 1:length(New_S3)
    New_fullV(end+1:end+length(New_S3{i}(:,37))) = New_S3{i}(:,37);
    Old_fullV(end+1:end+length(Old_S3{i}(:,37))) = Old_S3{i}(:,37);
end