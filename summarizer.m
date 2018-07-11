load Final_Models_Out

WT_time{1} = [0;T1{1}+200];
TG_time{1} = [0;T2{1}+200];
TGIso_time{1} = [0;T3{1}+200];

WT_V{1} = [S1{1}(1,37);S1{1}(:,37)];
TG_V{1} = [S2{1}(1,37);S2{1}(:,37)];
TGIso_V{1} = [S3{1}(1,37);S3{1}(:,37)];

WT_Ca{1} = [S1{1}(1,33);S1{1}(:,33)];
TG_Ca{1} = [S2{1}(1,33);S2{1}(:,33)];
TGIso_Ca{1} = [S3{1}(1,33);S3{1}(:,33)];

WT_INa{1} = [loginfo1.currents{1}(1,1);loginfo1.currents{1}(:,1)];
TG_INa{1} = [loginfo2.currents{1}(1,1);loginfo2.currents{1}(:,1)];
TGIso_INa{1} = [loginfo3.currents{1}(1,1);loginfo3.currents{1}(:,1)];

WT_ICa{1} = [loginfo1.currents{1}(1,11);loginfo1.currents{1}(:,11)];
TG_ICa{1} = [loginfo2.currents{1}(1,11);loginfo2.currents{1}(:,11)];
TGIso_ICa{1} = [loginfo3.currents{1}(1,11);loginfo3.currents{1}(:,11)];

WT_INaCa{1} = [loginfo1.currents{1}(1,14);loginfo1.currents{1}(:,14)];
TG_INaCa{1} = [loginfo2.currents{1}(1,14);loginfo2.currents{1}(:,14)];
TGIso_INaCa{1} = [loginfo3.currents{1}(1,14);loginfo3.currents{1}(:,14)];


for i = 2:50
    WT_time{i} = [0;T1{i}+200];
    TG_time{i} = [0;T2{i}+200];
    TGIso_time{i} = [0;T3{i}+200];

    WT_V{i} = [S1{i-1}(end,37);S1{i}(:,37)];
    TG_V{i} = [S2{i-1}(end,37);S2{i}(:,37)];
    TGIso_V{i} = [S3{i-1}(end,37);S3{i}(:,37)];

    WT_Ca{i} = [S1{i-1}(end,33);S1{i}(:,33)];
    TG_Ca{i} = [S2{i-1}(end,33);S2{i}(:,33)];
    TGIso_Ca{i} = [S3{i-1}(end,33);S3{i}(:,33)];
    
    WT_INa{1} = [loginfo1.currents{i-1}(end,1);loginfo1.currents{i}(:,1)];
    TG_INa{1} = [loginfo2.currents{i-1}(end,1);loginfo2.currents{i}(:,1)];
    TGIso_INa{1} = [loginfo3.currents{i-1}(end,1);loginfo3.currents{i}(:,1)];

    WT_ICa{1} = [loginfo1.currents{i-1}(end,11);loginfo1.currents{i}(:,11)];
    TG_ICa{1} = [loginfo2.currents{i-1}(end,11);loginfo2.currents{i}(:,11)];
    TGIso_ICa{1} = [loginfo3.currents{i-1}(end,11);loginfo3.currents{i}(:,11)];

    WT_INaCa{1} = [loginfo1.currents{i-1}(end,14);loginfo1.currents{i}(:,14)];
    TG_INaCa{1} = [loginfo2.currents{i-1}(end,14);loginfo2.currents{i}(:,14)];
    TGIso_INaCa{1} = [loginfo3.currents{i-1}(end,14);loginfo3.currents{i}(:,14)];   
end

subplot(2,1,1), plot(WT_time{50},WT_V{50},'b')
hold on 
plot(TG_time{50},TG_V{50},'g')
plot(TGIso_time{17},TGIso_V{17},'r')

subplot(2,1,2), plot(WT_time{50},WT_Ca{50},'b')
hold on 
plot(TG_time{50},TG_Ca{50},'g')
plot(TGIso_time{10},TGIso_Ca{10},'r')
pause
close

beats = [2,4,6,8,10,11,14,16];
ColOrd = linspecer(length(beats),'sequential');
[m,n] = size(ColOrd);
for i = 1:length(beats)
    subplot(2,1,1), plot(TGIso_time{beats(i)},TGIso_V{beats(i)},'Color',ColOrd{i})
    hold on 
    subplot(2,1,2), plot(TGIso_time{beats(i)},TGIso_Ca{beats(i)},'Color',ColOrd{i})
    hold on
    pause
end
colorbar

clearvars -except     WT_time TG_time...
    TGIso_time WT_V TG_V TGIso_V WT_Ca TG_Ca TGIso_Ca WT_INa TG_INa TGIso_INa...
    WT_ICa TG_ICa TGIso_ICa WT_INaCa TG_INaCa TGIso_INaCa

figure
plot(TGIso_time{10},TGIso_V{10},'k')
hold on
Final_Models_ablate_SR
Ablate_SR_time{1} = [0;T3{1}+200];
Ablate_SR{1} = [S3{1}(1,37);S3{1}(:,37)];
for i = 2:length(S3)
    Ablate_SR_time{i} = [0;T3{i}+200];
    Ablate_SR{i} = [S3{i-1}(end,37);S3{i}(:,37)];
end
plot(Ablate_SR_time{2},Ablate_SR{2},'r')

clearvars -except     WT_time TG_time...
    TGIso_time WT_V TG_V TGIso_V WT_Ca TG_Ca TGIso_Ca WT_INa TG_INa TGIso_INa...
    WT_ICa TG_ICa TGIso_ICa WT_INaCa TG_INaCa TGIso_INaCa Ablate_SR_time Ablate_SR...
    
Final_Models_ablate_INaCa
Ablate_INaCa_time{1} = [0;T3{1}+200];
Ablate_INaCa{1} = [S3{1}(1,37);S3{1}(:,37)];
for i = 2:length(S3)
    Ablate_INaCa_time{i} = [0;T3{i}+200];
    Ablate_INaCa{i} = [S3{i-1}(end,37);S3{i}(:,37)];
end
plot(Ablate_INaCa_time{2},Ablate_INaCa{2},'b')

clearvars -except     WT_time TG_time...
    TGIso_time WT_V TG_V TGIso_V WT_Ca TG_Ca TGIso_Ca WT_INa TG_INa TGIso_INa...
    WT_ICa TG_ICa TGIso_ICa WT_INaCa TG_INaCa TGIso_INaCa Ablate_SR_time Ablate_SR_V...
    Ablate_INaCa_time Ablate_INaCa

Final_Models_ablate_ICa
Ablate_ICa_time{1} = [0;T3{1}+200];
Ablate_ICa{1} = [S3{1}(1,37);S3{1}(:,37)];
for i = 2:length(S3)
    Ablate_ICa_time{i} = [0;T3{i}+200];
    Ablate_ICa{i} = [S3{i-1}(end,37);S3{i}(:,37)];
end
plot(Ablate_ICa_time{2},Ablate_ICa{2},'g')

clearvars -except     WT_time TG_time...
    TGIso_time WT_V TG_V TGIso_V WT_Ca TG_Ca TGIso_Ca WT_INa TG_INa TGIso_INa...
    WT_ICa TG_ICa TGIso_ICa WT_INaCa TG_INaCa TGIso_INaCa Ablate_SR_time Ablate_SR...
    Ablate_INaCa_time Ablate_INaCa
    