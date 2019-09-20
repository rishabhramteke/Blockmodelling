function mySpyPlot(A,mrksize,Z1,Z2,eta,Colors1,Colors2,PlotData)
% function to plot an adjacency matrix given the clustering by the IRM
% model
%
% A         I x J adjacency matrix
% mrksize   size of dots in spy-plot (default: 500/max(size(A{1},1)))
% Z1        clustering assignment matrix for mode 1
% Z2        clustering assignment matrix for mode 2
% eta       matrix of cluster relations
% Colors1   colors used for Z1
% Colors2   colors used for Z2
% PlotData  Boolean specifying wether to plot the links in the graph (default: true)

if ~iscell(A)
   B=A;
   clear A;
   A{1}=B;
end
if nargin<2 || isempty(mrksize)
    mrksize=500/max(size(A{1}));
end
if nargin<3 
    Z1=[];
end
if nargin<4 || isempty(Z2)
    Z2=Z1;
end
if nargin<5   
   eta=[];
end
if nargin<6 || isempty(Colors1)
   C1=colormap(lines);
   while size(Z1,1)>size(C1,1)
       C1=[C1; C1];
   end
   Colors1=Z1'*C1(1:size(Z1,1),:); 
end
if nargin<7 || isempty(Colors2)
   C2=colormap(lines);
   while size(Z2,1)>size(C2,1)
       C2=[C2; C2];
   end
   Colors2=Z2'*C2(1:size(Z2,1),:); 
end
if nargin<8
    PlotData = true;
end

nn=length(A);
for n=1:length(A)
    hold on;
    if ~isempty(eta)
        colormap(gray);
        IntensityImage = flipud(Z1'*log10(eta(:,:,n))*Z2);
        IntensityImage = (IntensityImage-min(IntensityImage(:))/range(IntensityImage(:)));
        IntensityImage(isinf(IntensityImage)) = 1/2;
        image(size(colormap,1)*(1-IntensityImage/2));
    end
    if nn>1
       subplot(ceil(sqrt(nn)),ceil(sqrt(nn)),n);
       title(['A\{' num2str(n) '\}'],'FontWeight','bold')
    end
    [I,J,val]=find(A{n});
    linecol = [.5 .5 1];
    hold on;
    [N1 N2]=size(A{n});
    if ~isempty(Z1)
        
        sZ1=cumsum(sum(Z1,2));
        sZ2=cumsum(sum(Z2,2));
        for k=1:length(sZ1)-1            
           plot([0.5 sZ2(end)+0.5], N1+1-[sZ1(k)+.5 sZ1(k)+.5],'-','LineWidth',0.5,'Color',linecol); 
        end
        for k=1:length(sZ2)-1            
           plot([sZ2(k)+.5 sZ2(k)+.5],N1+1-[0.5 sZ1(end)+.5],'-','LineWidth',0.5,'Color',linecol);
        end
    end
    
    if PlotData
        plot(J,N1+1-I,'.k','MarkerSize',mrksize);
    end
       
    plot([0.5 N2+0.5], [0.5 0.5],'-k','LineWidth',2);
    plot([0.5 0.5], [0.5 N1+0.5],'-k','LineWidth',2); 
    plot([N2+0.5 N2+0.5], [0.5 N1+0.5],'-k','LineWidth',2); 
    plot([0.5 N2+0.5], [N1+0.5 N1+0.5],'-k','LineWidth',2); 
    
    if nargin>2
        for v=1:size(Z1,2)
            plot([-0.5 -0.5],N1+1-[v-0.5 v+0.5],'-','color',Colors1(v,:),'LineWidth',4);             
        end    
        for v=1:size(Z2,2)
            plot([v-0.5 v+0.5],[N1+1.5 N1+1.5],'-','color',Colors2(v,:),'LineWidth',4);
        end
    end
    axis equal;
    axis([-1 N2+5 -1 N1+5]);
    axis off;    
    set(gca,'XTick',[])
    set(gca,'XTickLabel',[])
    set(gca,'YTick',[])
    set(gca,'YTickLabel',[])
end