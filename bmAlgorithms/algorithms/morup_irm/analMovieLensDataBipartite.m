clear all;
close all;
load MovielensData

%% Analyze Graph as Bipartite Network
[W,class]=createValidationData(A,2.5,'Bipartite');
noc=[10 10];
opts.maxiter=20;
opts.nsampleiter=10;
opts.type='Binary';
[L,cpu_time,Zuser,Zmovie,eta,sample,West]=IRMBipartite(A,W,noc,opts);

%% Evaluate the predictive performance of the model
[AUC,TPR,FPR]=calcAUC(West,A);

%% Display logP, i.e. the value of the join log-posterior log P(A,Z|\alpha,beta^+,beta^-)
subplot(1,2,1);
plot(L);
xlabel('Iteration','FontWeight','bold')
ylabel('logP','FontWeight','bold');
% Display the ROC for the prediction of links according to the model
subplot(1,2,2);
hold on;
plot(FPR,TPR,'LineWidth',2);
plot([0 1],[0 1],':k','LineWidth',2);
xlabel('FPR','FontWeight','bold')
ylabel('TPR','FontWeight','bold');
title(['AUC score ' num2str(AUC)],'FontWeight','bold');

%% Sort the graph and display the sorted results
[Asorted,Zuser_sorted,Zmovie_sorted,eta_sorted,perm1,perm2]=sortGraphBipartite(A,full(Zuser),full(Zmovie),eta);
figure; 
subplot(2,6,1:3); 
mySpyPlot(A); 
title('Unsorted graph','FontWeight','Bold'); 
subplot(2,6,4:6); 
mySpyPlot(Asorted,[],Zuser_sorted,Zmovie_sorted,eta_sorted); 
title('Graph sorted according to extracted assignments','FontWeight','Bold');
subplot(2,6,7:8)
imagesc(-Zuser_sorted);
title('Estimated User Groups (Zuser_{sorted})','FontWeight','bold')
subplot(2,6,9:10)
imagesc(-eta_sorted); colorbar;
axis equal;
axis tight;
xlabel('Movie cluster nr','FontWeight','Bold')
ylabel('User cluster nr','FontWeight','Bold')
title('\eta_{sorted}','FontWeight','bold');
subplot(2,6,11:12)
imagesc(-Zmovie_sorted);
title('Estimated Movie Groups (Zmovie_{sorted})','FontWeight','bold')

%% Display movies in clusters
figure('Name','Examples of up to 20 movies from each extracted movie cluster');
Names=movie{2}(perm2);
for d=1:size(Zmovie_sorted,1);
   ind=find(Zmovie_sorted(d,:)==1);
   subplot(3,ceil(size(Zmovie_sorted,1)/3),d);
   for i=1:min([20,length(ind)])
       text(0,1-(i-1)/20,Names{ind(i)}(1:min([length(Names{ind(i)}),20])),'FontSize',6);    
   end          
   axis off;
   title(['Cluster ' num2str(d)],'FontWeight','bold')
end

%% Summarize the extracted movie groups in terms of genre informatioin
figure('Name','Distribution of movie genres in each extracted movie group')
m_genre=movie{6}(perm2,:);
for d=1:size(Zmovie_sorted,1);
   ind=find(Zmovie_sorted(d,:)==1);
   subplot(5,ceil(size(Zmovie_sorted,1)/5),d);
   N=sum(m_genre(ind(i),:),1);      
   bar(genre{2},N/sum(N));      
   axis([-0.5 18.5 0 1]);
   title(['Cluster ' num2str(d)],'FontWeight','bold')
end
disp('Interpretation of Genre:')
for k=1:length(genre{1})
    disp([num2str(genre{2}(k)) ': ' genre{1}{k}])
end

%% Summary of user groups
figure('Name','Distribution of gender and age for each extracted user group')
clear N;
for d=1:size(Zuser_sorted,1);
   ind=find(Zuser_sorted(d,:)==1);
   subplot(5,ceil(size(Zuser_sorted,1)/5),d);
   % Female and Male Gender
   N(1)=sum(strcmp(user{3}(perm1(ind)),'M'))/length(ind);
   N(2)=1-N(1);
   % Age groups <20 20-30 30-40 >=40
   N(3)=sum(user{2}(perm1(ind))<20)/length(ind);
   N(4)=sum(user{2}(perm1(ind))>=20 & user{2}(perm1(ind))<30)/length(ind);
   N(5)=sum(user{2}(perm1(ind))>=30 & user{2}(perm1(ind))<40)/length(ind);
   N(6)=sum(user{2}(perm1(ind))>=40)/length(ind);      
   hold on;
   h=bar(N);   
   plot([2.5 2.5],[0 1],'r-','LineWidth',2);
   axis([0.5 6.5 0 1]);
   title(['Cluster ' num2str(d)],'FontWeight','bold')
end
disp('Interpretation of Users:')
disp('1 :  Male')
disp(['2 :  Female'])
disp(['3 :  <20 years'])
disp(['4 :  20-29 years'])
disp(['5 :  30-39 years'])
disp(['6 :  >=40 years'])



