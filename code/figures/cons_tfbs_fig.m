%Plot cons scores and TFBS Ranking

%import cons scores
[cons_scores ~] = importdata('C:\Users\Nicholas\Documents\GitHub\synth_es2\data\analysis\ToFIMO\Results_12.01.16\SegmentedES2_CONSscores.csv');
cons_scores = cons_scores.data;
scores = cons_scores(:,2);

%import tfbs results from fimo
tfbs = csvread('C:\Users\Nicholas\Documents\GitHub\synth_es2\data\analysis\ToFIMO\Results_12.01.16\TFBSfiltered_plot.csv');

%find joint p value per location
tfbs_scores = ones(length(cons_scores),1);
for i = 1:size(tfbs,1)
    start = tfbs(i,1);
    stop = tfbs(i,2);
    for j = start:stop
        tfbs_scores(j) = min(tfbs_scores(j),tfbs(i,3));
    end
end
tfbs_scores = ones(length(tfbs_scores),1)-tfbs_scores-min(tfbs_scores(tfbs_scores > 0))*ones(length(tfbs_scores),1);
tfbs_scores = tfbs_scores./max(tfbs_scores);


%plot cons and tfbs together

% x = 1:length(tfbs_scores);
% plot(x,tfbs_scores,x,cons_scores(:,2),'black', 'LineWidth',2);
figure;
area(cons_scores(:,2));

% plot(tfbs_scores, 'LineWidth',2);
% legend({'pValue','Conservation'}, 'Location', 'southeast');
title('Conservation Scores');

ylabel('Normalized Score');

figure;
area(tfbs_scores, 'FaceColor', 'black');
title('TFBS pValues');
xlabel('Location (bp)');
ylabel('Normalized Score');


 axis([0 900 0 1]);
