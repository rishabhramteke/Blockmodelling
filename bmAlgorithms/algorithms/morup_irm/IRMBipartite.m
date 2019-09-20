function [L,cpu_time,Z1,Z2,eta,sample,West]=IRMBipartite(A,W,noc,opts)

% Non-parametric IRM of Bi-partite multi-graphs based on collapsed Gibbs sampling with split
% merge moves
%
% Usage: 
%  [L,cpu_time,Z1,Z2,eta,sample,West]=IRMBipartite(A,W,noc,opts)
%
% Input:
%   A       cell array of I x J  sparse adjacency matrices
%   W       cell array of I x J  sparse missing value indicator matrices 
%   noc     1 x 2  vector indicating the number of initial clusters in Z1 and Z2
%   opts.
%           maxiter     maximum number of iterations
%           Z1           noc x I matrix of community identities
%           Z2           noc x J matrix of community identities
%           dZstep      number of iteratiosn between each recorded sample
%           verbose     1: display iteration results, 0: no display
%           eta0        1x2 vector of pseudo link and non-link counts between the groups (default: [1 1])
%           type        'Binary' (default), 'Weighted'
% Output:
%   L           Log-likelihood function at each iteration
%   cpu_time    cpu-time cost for each iteration
%   Z1          Estimated row-clustering assigment matrix
%   Z2          Estimated column-clustering assigment matrix
%   sample      sampled parameters at each dZstep iterationCluster
%   eta         Average between group link densitites
%   West        Estimated link values for missing links and non-links
%               for type='Categorical' this corresponds to the probability
%               of generating the observed class.
 
if nargin<4
    opts=struct;
end

% If A is adjacency matrix put it in cell structure
if ~iscell(A)
    B=A;
    clear A;
    A{1}=B;
    clear B;
end
if ~iscell(W)
    B=W;
    clear W;
    W{1}=B;
    clear B;
end
Ngraphs=length(A);
nnzW=0;
% Make sure links treated as missing are removed from estimation
par.type=mgetopt(opts,'type','Binary');
par.categories=[];
[I,J]=size(A{1});
for k=1:Ngraphs
    W{k}=logical(W{k});
    V=W{k}.*A{k}+W{k};    
    West{k}=sparse(I,J);
    [Iw{k},Jw{k},vv{k}]=find(V);
    vv{k}=vv{k}-1;
    clear V;    
    if strcmp(par.type,'Binary')
        A{k}=logical(A{k}-A{k}.*W{k}); % Remove missing links
    else
         A{k}=A{k}-A{k}.*W{k}; % Remove missing links
    end
    par.categories=union(par.categories,full(unique(A{k})));
    nnzW=nnzW+nnz(W{k});
end
par.categories=setdiff(par.categories,0);

% Initialize Z1 and Z2
ind=ceil(noc(1)*rand(1,I));
Z1=sparse(ind,1:I,ones(1,I),noc(1),I);
Z1=mgetopt(opts,'Z1',Z1);
Z1=full(Z1);
Z1(sum(Z1,2)==0,:)=[];
ind=ceil(noc(2)*rand(1,J));
Z2=sparse(ind,1:J,ones(1,J),noc(2),J);
Z2=mgetopt(opts,'Z2',Z2);
Z2=full(Z2);
Z2(sum(Z2,2)==0,:)=[];
% Zet remaining parameters
alpha1=mgetopt(opts,'alpha1',log(I));
alpha2=mgetopt(opts,'alpha2',log(J));
maxiter=mgetopt(opts,'maxiter',50);
verbose=mgetopt(opts,'verbose',1);
switch par.type
    case 'Binary'
        eta0=mgetopt(opts,'eta0',[1 1]); % pseudo counts of links and non-links between clusters
    case 'Weighted'
        eta0=mgetopt(opts,'eta0',[1 1e-6]); % pseudo counts of links and non-links between clusters
end
sample_step=mgetopt(opts,'sample_step',25); % Zteps between samples
nsampleiter=mgetopt(opts,'nsampleiter',25);

plotfcn = mgetopt(opts, 'plotfcn', []);

% Zet algorithm variables
sample=struct;
L=zeros(1,maxiter);
cpu_time=L;

sstep=0;
westiter=0;
Q=-inf;
Qbest=-inf;
iter=0;

if verbose % Display algorithm    
    disp(['Non-parametric bi-partite clustering based on the IRM model for ' par.type ' graphs'])
    dheader = sprintf('%12s | %12s | %12s | %12s | %12s | %12s | %12s','Iteration','logP','dlogP/|logP|','noc1','noc2','time');
    dline = sprintf('-------------+--------------+--------------+--------------+--------------+--------------');
    disp(dline);
    disp(dheader);
    disp(dline);
end


% Main Loop
while iter<maxiter
   
    iter=iter+1;
    tic;
    Qold=Q;
    
    % Gibbs sampling of Z1 
    [Z1,logP_A,logP_Z1]=Gibbs_sample_ZIRM(Z1,Z2,A,W,eta0,alpha1,randperm(I),1,par);    
    for t=1
        [Z1,logP_A,logP_Z1]=split_merge_sample_Z(Z1,Z2,A,W,eta0,alpha1,logP_A,logP_Z1,1,par);
    end    
    noc(1)=size(Z1,1);
    
    % Gibbs sampling of Z2 
    [Z2,logP_A,logP_Z2]=Gibbs_sample_ZIRM(Z2,Z1,A,W,eta0,alpha2,randperm(J),2,par);        
    for t=1
        [Z2,logP_A,logP_Z2]=split_merge_sample_Z(Z2,Z1,A,W,eta0,alpha2,logP_A,logP_Z2,2,par);    
    end
    noc(2)=size(Z2,1);
    
    % Evaluate result
    Q=logP_A+logP_Z1+logP_Z2;
    dQ=Q-Qold;    
    L(iter)=Q;
    t_iter=toc;
    cpu_time(iter)=t_iter;    
   
    % Display iteration
    if rem(iter,1)==0 && verbose        
        disp(sprintf('%12.0f | %12.4e | %12.4e | %12.0f | %12.0f |%12.4f ',iter,Q,dQ/abs(Q),noc(1),noc(2),t_iter));
    end
    eta=calculateRho(Z2,Z1,A,W,eta0,par);     
    
    % Ztore sample
    if mod(iter,sample_step)==0         
        sstep=sstep+1;
        sample.iteration(sstep)=iter;
        sample.Z1{sstep}=Z1;        
        sample.Z2{sstep}=Z2;        
        sample.eta{sstep}=eta;
    end        
    if Q>Qbest
       sample.MAP.Z1=Z1; 
       sample.MAP.Z2=Z2;
       sample.MAP.iter=iter;
       sample.MAP.L=Q;
       sample.MAP.eta=eta;
       Qbest=Q;
    end
        
        
    % Estimate missing link probability of sample
    if iter>maxiter-nsampleiter && nnzW>0
        westiter=westiter+1;        
        if iter==maxiter-nsampleiter+1
            disp(['Initiating estimation of missing links for the last ' num2str(nsampleiter) ' iterations']);   
        end      
        step=10000;
        for k=1:Ngraphs
            val=zeros(1,length(Iw{k}));            
            switch par.type
                case {'Binary','Weighted'}
                    for kk=1:ceil((length(Iw{k})/step))
                        ind=(kk-1)*step+1:min([kk*step, length(Iw{k})]);   
                        val(ind)=sum(Z1(:,Iw{k}(ind)).*(eta(:,:,k)*Z2(:,Jw{k}(ind))))+eps;
                    end    
                case 'Categorical'
                    for c=0:length(par.categories)
                        indc=(vv{k}==c);
                        for kk=1:ceil((length(indc)/step))
                            ind=(kk-1)*step+1:min([kk*step, length(indc)]);
                            if c==0
                                 val(indc(ind))=sum(Z1(:,Iw{k}(indc(ind))).*((1-sum(eta(:,:,k,:),4))*Z2(:,Jw{k}(indc(ind)))))+eps;
                            else
                                 val(indc(ind))=sum(Z1(:,Iw{k}(indc(ind))).*(eta(:,:,k,c)*Z2(:,Jw{k}(indc(ind)))))+eps;
                            end
                        end                                        
                    end
            end
            West{k}=West{k}+sparse(Iw{k},Jw{k},val,size(A{k},1),size(A{k},2));
        end
    end
    
    % plot result
    if isa(plotfcn, 'function_handle'), plotfcn(); end;
end

% Average link predictions
if nnzW>0
    for k=1:Ngraphs
        West{k}=West{k}/westiter;
    end
end
[v,ind1]=sort(sum(Z1,2),'descend');
[v,ind2]=sort(sum(Z2,2),'descend');
Z1=Z1(ind1,:);
eta=eta(ind1,:,:,:);
Z2=Z2(ind2,:);
eta=eta(:,ind2,:,:);

Z1=sparse(Z1);
Z2=sparse(Z2);

% Display final iteration
if verbose   
  disp('Result of final iteration');
  disp(sprintf('%12.0f | %12.4e | %12.4e | %12.0f | %12.0f |%12.4f ',iter,Q,dQ/abs(Q),noc(1),noc(2),t_iter));
end

% -------------------------------------------------------------------------
% Parser for optional arguments
function var = mgetopt(opts, varname, default, varargin)
if isfield(opts, varname)
    var = getfield(opts, varname); 
else
    var = default;
end
for narg = 1:2:nargin-4
    cmd = varargin{narg};
    arg = varargin{narg+1};
    switch cmd
        case 'instrset',
            if ~any(strcmp(arg, var))
                fprintf(['Wrong argument %s = ''%s'' - ', ...
                    'Using default : %s = ''%s''\n'], ...
                    varname, var, varname, default);
                var = default;
            end
        otherwise,
            error('Wrong option: %s.', cmd);
    end
end

%-------------------------------------------------------------------------
function eta=calculateRho(Z2,Z1,A,W,eta0,par)
        noc1=size(Z1,1);
        noc2=size(Z2,1);
        sumZ1=sum(Z1,2);
        sumZ2=sum(Z2,2);
        Ngraphs=length(A);
        if strcmp(par.type,'Categorical')
            eta=zeros(noc1,noc2,Ngraphs,length(par.categories));
        else
            eta=zeros(noc1,noc2,Ngraphs);
        end
        for k=1:Ngraphs            
            switch par.type
                case 'Binary'
                    ZAZt=Z1*A{k}*Z2';
                    ZWZt=Z1*W{k}*Z2';            
                    ZZ=sumZ1*sumZ2'-ZAZt-ZWZt;                
                    n_link=ZAZt+eta0(1);                        
                    n_nonlink=ZZ+eta0(2);   
                    eta(:,:,k)=n_link./(n_link+n_nonlink);        
                case 'Categorical'
                    ZAZt=zeros(noc1,noc2,length(par.categories));
                    for c=1:length(par.categories)
                        ZAZt(:,:,c)=Z1*(A{k}==par.categories(c))*Z2';
                    end
                    ZWZt=Z1*W{k}*Z2';            
                    ZZ=sumZ1*sumZ2'-sum(ZAZt,3)-ZWZt;                
                    n_link=ZAZt+eta0(1);                        
                    n_nonlink=ZZ+eta0(2);   
                    eta(:,:,k,:)=n_link./repmat((sum(n_link,3)+n_nonlink),[1 1 length(par.categories)]);        
                case 'Weighted'
                    ZAZt=Z1*A{k}*Z2';
                    ZWZt=Z1*W{k}*Z2';            
                    ZZ=sumZ1*sumZ2'-ZWZt;                
                    n_link=ZAZt+eta0(1);                        
                    n_nonlink=ZZ+eta0(2);   
                    eta(:,:,k)=n_link./(n_nonlink);        
            end
        end

% -------------------------------------------------------------------------  
function [Z2,logP_A,logP_Z2]=split_merge_sample_Z(Z2,Z1,A,W,eta0,alpha,logP_A,logP_Z2,mode,par)
    %[logP_A_t,logP_Z2_t]=evalLikelihood(Z2,Z1,A,W,eta0,alpha,mode);                         
    %logP_A
    %logP_A_t
    %logP_Z2
    %logP_Z2_t
    

    noc2=size(Z2,1);
    if mode==2
        J=size(A{1},2);    
    else
        J=size(A{1},1);    
    end
    
    % step 1 select two observations i and j        
    ind1=ceil(J*rand);        
    ind2=ceil((J-1)*rand);
    if ind1<=ind2
       ind2=ind2+1;
    end
    clust1=find(Z2(:,ind1));
    clust2=find(Z2(:,ind2));

    if clust1==clust2 % Split   
        setZ=find(sum(Z2([clust1 clust2],:)));    
        setZ=setdiff(setZ,[ind1,ind2]);
        n_setZ=length(setZ);
        Z2_t=Z2;
        Z2_t(clust1,:)=0;        
        comp=[clust1 noc2+1];               
        Z2_t(comp(1),ind1)=1;
        Z2_t(comp(2),ind2)=1;

        % Reassign by restricted gibbs sampling        
        if n_setZ>0
            for rep=1:3
                [Z2_t,logP_A_t,logP_Z2_t,logQ_trans,comp]=Gibbs_sample_ZIRM(Z2_t,Z1,A,W,eta0,alpha,setZ(randperm(n_setZ)),mode,par,comp);                     
            end     
        else
           logQ_trans=0;
           [logP_A_t,logP_Z2_t]=evalLikelihood(Z2_t,Z1,A,W,eta0,alpha,mode,par);                 
        end
        
        % Calculate Metropolis-Hastings ratio
        a_split=rand<exp(logP_A_t+logP_Z2_t-logP_A-logP_Z2-logQ_trans);         
        if a_split
           disp(['Splitting cluster ' num2str(clust1) ' in mode ' num2str(mode)])
           logP_A=logP_A_t;
           logP_Z2=logP_Z2_t;
           Z2=Z2_t;
        end
    else % Merge                                     
        Z2_t=Z2;
        Z2_t(clust1,:)=Z2_t(clust1,:)+Z2_t(clust2,:);
        setZ=find(Z2_t(clust1,:));           
        Z2_t(clust2,:)=[];        
        if clust2<clust1
            clust1_t=clust1-1;
        else 
            clust1_t=clust1;
        end
        noc2_t=noc2-1;

        % calculate likelihood of merged cluster       
        [logP_A_t,logP_Z2_t]=evalLikelihood(Z2_t,Z1,A,W,eta0,alpha,mode,par);                

        % Zplit the merged cluster and calculate transition probabilties                        
        setZ=setdiff(setZ,[ind1,ind2]);
        n_setZ=length(setZ);
        Z2_tt=Z2_t;
        Z2_tt(clust1_t,:)=0;        
        comp=[clust1_t noc2_t+1];               
        Z2_tt(comp(1),ind1)=1;
        Z2_tt(comp(2),ind2)=1;                
        
        % Reassign by restricted gibbs sampling
        if n_setZ>0
            for rep=1:2        
                [Z2_tt,logP_A_tt,logP_Z2_tt,logQ_trans,comp]=Gibbs_sample_ZIRM(Z2_tt,Z1,A,W,eta0,alpha,setZ(randperm(n_setZ)),mode,par,comp);                
            end
            Force=[1 2]*Z2([clust1 clust2],:);        
            [Z2_tt,logP_A_tt,logP_Z2_tt,logQ_trans]=Gibbs_sample_ZIRM(Z2_tt,Z1,A,W,eta0,alpha,setZ(randperm(n_setZ)),mode,par,comp,Force);                        
        else
            logQ_trans=0;                   
        end
        a_merge=rand<exp(logP_A_t+logP_Z2_t-logP_A-logP_Z2+logQ_trans); 
        
        if a_merge
          disp(['Merging cluster ' num2str(clust1) ' with cluster ' num2str(clust2) ' in mode ' num2str(mode)])
          logP_A=logP_A_t;
          logP_Z2=logP_Z2_t;
          Z2=Z2_t;          
        end
    end

% -------------------------------------------------------------------------  
function [logP_A,logP_Z2,ZAZt,ZWZt,sumZ2]=evalLikelihood(Z2,Z1,A,W,eta0,alpha,mode,par)

    switch par.type
        case 'Binary'
            my_func=@betaln;
            const=my_func(eta0(1),eta0(2));   
        case 'Weighted'
            my_func=@poissonln;
            const=my_func(eta0(1),eta0(2));           
    end
    if mode==2
        J=size(A{1},2);
    else
        J=size(A{1},1);
    end
    noc1=size(Z1,1);
    noc2=size(Z2,1);
    Ngraphs=length(A);
    ZAZt=zeros(noc1,noc2,Ngraphs);
    ZWZt=ZAZt;
            
    for k=1:length(A)
        if mode==2
            ZAZt(:,:,k)=Z1*A{k}*Z2';
            ZWZt(:,:,k)=Z1*W{k}*Z2';
        else
            ZAZt(:,:,k)=Z1*A{k}'*Z2';
            ZWZt(:,:,k)=Z1*W{k}'*Z2';
        end    
    end
    sumZ1=sum(Z1,2);
    sumZ2=sum(Z2,2);
    switch par.type
        case 'Binary'
            ZZ=repmat(sumZ1*sumZ2',[1 1 Ngraphs])-ZAZt-ZWZt;                
        case 'Weighted'
            ZZ=repmat(sumZ1*sumZ2',[1 1 Ngraphs])-ZWZt;                        
    end
    n_link=ZAZt+eta0(1);    
    n_nonlink=ZZ+eta0(2);         
    logP_A=sum(sum(sum(my_func(n_link,n_nonlink))))-Ngraphs*noc1*noc2*const;
    logP_Z2=noc2*log(alpha)+sum(gammaln(full(sumZ2)))-gammaln(J+alpha)+gammaln(alpha);      
    
% -------------------------------------------------------------------------
function [Z2,logP_A,logP_Z2,logQ_trans,comp]=Gibbs_sample_ZIRM(Z2,Z1,A,W,eta0,alpha,JJ,mode,par,comp,Force)    
    
    if nargin<11
        Force=[];
    end
    if nargin<10
        comp=[];
    end
    logQ_trans=0;

    switch par.type
        case 'Binary'
            my_func=@betaln; 
            const=my_func(eta0(1),eta0(2));                   
        case 'Weighted' 
            my_func=@poissonln;
            const=my_func(eta0(1),eta0(2));                           
    end
    if mode==2
        [I,J]=size(A{1});
    else
        [J,I]=size(A{1});
    end
    Ngraphs=length(A);    
        
    t=0;   
    sumZ1=sum(Z1,2);
    noc1=length(sumZ1);    
    sumZ2=sum(Z2,2);
    noc2=length(sumZ2);        
    
    ZA=zeros(noc1,J,Ngraphs);
    ZW=zeros(noc1,J,Ngraphs);
    n_link=zeros(noc1,noc2,Ngraphs);
    n_nonlink=n_link;
    
    for k=1:Ngraphs     
        if mode==2
            ZA(:,:,k)=Z1*A{k};
            ZW(:,:,k)=Z1*W{k};
        else
            ZA(:,:,k)=Z1*A{k}';
            ZW(:,:,k)=Z1*W{k}';
        end
        ZAZt=ZA(:,:,k)*Z2';
        ZWZt=ZW(:,:,k)*Z2';                    
        n_link(:,:,k,:)=ZAZt+eta0(1);
        switch par.type
            case 'Binary'
                n_nonlink(:,:,k)=sumZ1*sumZ2'-sum(ZAZt,3)-ZWZt+eta0(2);        
            case 'Weighted'
                n_nonlink(:,:,k)=sumZ1*sumZ2'-ZWZt+eta0(2);        
        end
    end    
    beta_eval=my_func(n_link,n_nonlink);
    sumZ1Ngraphs=repmat(sumZ1,[1,1,Ngraphs]);
    
    for k=JJ           
        t=t+1;
        if mod(t,5000)==0
            disp(['sampling ' num2str(t) ' out of ' num2str(J) ' nodes']);
        end
        
        % Remove effect of Z2(:,k)        
        sumZ2=sumZ2-Z2(:,k);   
        Z1AZ2k=ZA(:,k,:,:);
        Z1WZ2k=ZW(:,k,:);  
        switch par.type
            case 'Binary'
                nZ1AZ2k=sumZ1Ngraphs-sum(Z1AZ2k,4)-Z1WZ2k;
            case 'Weighted'
                nZ1AZ2k=sumZ1Ngraphs-Z1WZ2k;
        end
        d=find(Z2(:,k));        
        if ~isempty(d)
            n_link(:,d,:,:)=n_link(:,d,:,:)-Z1AZ2k;                                            
            n_nonlink(:,d,:)=n_nonlink(:,d,:)-nZ1AZ2k;                                                       
            Z2(:,k)=0;               
        end
        
        if isempty(comp)
            if sumZ2(d)==0 % Remove singleton cluster            
                v=1:noc2;
                v(d)=[];
                d=[];
                noc2=noc2-1;                               
                P=sparse(1:noc2,v,ones(1,noc2),noc2,noc2+1);                        
                Z2=P*Z2;                        
                sumZ2=sumZ2(v,1);            
                n_link=n_link(:,v,:,:);            
                n_nonlink=n_nonlink(:,v,:);            
                beta_eval=beta_eval(:,v,:);            
            end                                       

            % Calculate probability for existing communties as well as proposal cluster                                                                 
            beta_eval(:,d,:)=my_func(n_link(:,d,:,:),n_nonlink(:,d,:)); % removed the constant -my_func(Ap,An)))                                           
            sum_beta_eval=sum(sum(sum(beta_eval)));
            e=ones(1,noc2);
            sum_beta_eval_d=sum_beta_eval-[sum(sum(beta_eval,1),3) 0];
            TTT1=cat(2,n_link+Z1AZ2k(:,e,:,:),Z1AZ2k+eta0(1));
            TTT2=cat(2,n_nonlink+nZ1AZ2k(:,e,:),nZ1AZ2k+eta0(2));
            beta_eval_d=my_func(TTT1,TTT2);                                        
            logQ=sum_beta_eval_d'+[sum(sum(beta_eval_d(:,1:noc2,:),1),3) (sum(sum(beta_eval_d(:,noc2+1,:),1),3)-Ngraphs*noc1*const)]'; % removed the constant -my_func(Ap,An)))                                             

            % Zample from posterior                        
            QQ=exp(logQ-max(logQ));
            weight=[sumZ2; alpha];
            QQ=weight.*QQ;        
            ind=find(rand<cumsum(QQ/sum(QQ)),1,'first');     
            Z2(ind,k)=1;   
            if ind>noc2            
                noc2=noc2+1;
                sumZ2(noc2,1)=0;
                n_link(:,noc2,:,:)=eta0(1);                                
                n_nonlink(:,noc2,:)=eta0(2);                     
                beta_eval(:,noc2,:)=0;                    
            end
            beta_eval_d=beta_eval_d(:,ind,:);            
        else            
            % Calculate probability for existing communties as well as proposal cluster                                                                  
            beta_eval(:,d,:)=my_func(n_link(:,d,:,:),n_nonlink(:,d,:)); % removed the constant -my_func(Ap,An)))                                           
            sum_beta_eval=sum(sum(sum(beta_eval)));                        
            e=ones(1,2);
            sum_beta_eval_d=sum_beta_eval-sum(sum(beta_eval(:,comp,:),1),3);
            beta_eval_d=my_func(n_link(:,comp,:,:)+Z1AZ2k(:,e,:,:),n_nonlink(:,comp,:)+nZ1AZ2k(:,e,:));                                        
            logQ=sum_beta_eval_d'+sum(sum(beta_eval_d,1),3)'; % removed the constant -my_func(Ap,An)))                                             
           
            % Zample from posterior                        
            QQ=exp(logQ-max(logQ));
            weight=sumZ2(comp);
            QQ=weight.*QQ;
            QQ=QQ/sum(QQ);            
            if isempty(Force)
                ind=find(rand<cumsum(QQ),1,'first');
            else 
                ind=Force(k);
            end
            logQ_trans=logQ_trans+log(QQ(ind)+eps);
            Z2(comp(ind),k)=1;  
            beta_eval_d=beta_eval_d(:,ind,:);            
            ind=comp(ind);
        end
                        
        % Re-enter effect of new s_k        
        sumZ2=sumZ2+Z2(:,k);
        n_link(:,ind,:,:)=n_link(:,ind,:,:)+Z1AZ2k;                                
        n_nonlink(:,ind,:)=n_nonlink(:,ind,:)+nZ1AZ2k;        
        beta_eval(:,ind,:)=beta_eval_d;
                
        % Remove empty clusters        
        if ~all(sumZ2)
            d=find(sumZ2==0);
            if ~isempty(comp)
                ind_d=find(d<comp);
                comp(ind_d)=comp(ind_d)-1;
            end
            v=1:noc2;
            v(d)=[];
            noc2=noc2-length(d);                               
            P=sparse(1:noc2,v,ones(1,noc2),noc2,noc2+1);                        
            Z2=P*Z2;                        
            sumZ2=sumZ2(v,1);            
            n_link=n_link(:,v,:,:);
            n_nonlink=n_nonlink(:,v,:);
            beta_eval=beta_eval(:,v,:);            
        end              
    end      
    % Calculate Likelihood for sampled solution     
    logP_Z2=noc2*log(alpha)+sum(gammaln(full(sumZ2)))-gammaln(J+alpha)+gammaln(alpha);
    logP_A=sum(sum(sum(beta_eval)))-Ngraphs*noc1*noc2*const;     
    
    % ---------------------
    function C=poissonln(A,B)
        C = gammaln(A)-A.*log(B);
        
        