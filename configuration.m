function [data,pos,n_config]=configuration(type,l,max_ze,n_max,plotit)
%% CONFIGURATION compute the position of different configuration.
% 
% 
% * type        defining the type of configuration 'wenner', 'pole-pole',
%               'dipole-dipole',...
% * l           number of electrode [no unit]
% * plotit      Bolean (1,0) decide for ploting or not
% 


%%
% * Defining for the selected type
switch type
    case 'wenner'
        A_fx    = @(u,n,a) u;
        M_fx    = @(u,n,a) u+a;
        N_fx    = @(u,n,a) u+a*2;
        B_fx    = @(u,n,a) u+a*3;
        good    = @(u,n,a,l) u+(a+1)*3<=l && n<=1;
        ze_fx   = @(u,n,a) 0.173*a*3;
        xa_fx   = @(u,n,a) u+a*1.5;
        k_fx    = @(n,a) 2*pi*a;
    case 'pole-pole'
        A_fx    = @(u,n,a) u;
        M_fx    = @(u,n,a) u+a;
        N_fx    = @(u,n,a) NaN;
        B_fx    = @(u,n,a) NaN;
        good    = @(u,n,a,l) u+(a+1)<=l && n<=1;
        ze_fx   = @(u,n,a) 0.35*a;
        xa_fx   = @(u,n,a) u+a*0.5;
        k_fx    = @(n,a) 2*pi*a;
    case 'wenner-schlumberger'
        zec=[.173 .186 .189 .190];
        A_fx    = @(u,n,a) u;
        M_fx    = @(u,n,a) u+a*n;
        N_fx    = @(u,n,a) u+a*n+a;
        B_fx    = @(u,n,a) u+a*n+a+a*n;
        good    = @(u,n,a,l) (u+(a+1)*n+(a+1)+(a+1)*n)<=l;
        ze_fx   = @(u,n,a) zec(min(n,length(zec)))*(a*n+a+a*n);
        xa_fx   = @(u,n,a) u+n*a+a*0.5;
        k_fx    = @(n,a) pi*n*(n+1)*a;
    case 'dipole-dipole'
        zec=[.139 .174 .192 .203 .211 .216 .22 .224 .255];
        B_fx    = @(u,n,a) u;
        A_fx    = @(u,n,a) u+a;
        M_fx    = @(u,n,a) u+a+n*a;
        N_fx    = @(u,n,a) u+a+a*n+a;
        good    = @(u,n,a,l) (u+(a+1)+(a+1)*n+(a+1))<=l;
        ze_fx   = @(u,n,a) zec(min(n,length(zec)))*(a+a*n+a);
        xa_fx   = @(u,n,a) u+a+n*a*0.5;
        k_fx    = @(n,a)   pi*n*(n+1)*(n+2)*a;
        
end


%%
% Computing all possible configuration
i=0;
u=0;
a=0;
good_u=true;
while good_u
    u=u+1;
    n=0;
    if ~good(u,n,a,l)
        good_u=false;
    else
        good_n=true;
    end
    while good_n
        n=n+1;
        a=0;
        if ~good(u,n,a,l)
            good_n=false;
        else
            good_a=true;
        end
        while good_a
            a=a+1;
            i=i+1;
            data(i,:)   = [A_fx(u,n,a) B_fx(u,n,a) M_fx(u,n,a) N_fx(u,n,a)];
            pos(i,:)    = [xa_fx(u,n,a) ze_fx(u,n,a)];
            k(i)        = k_fx(n,a);
            assert(any(data(i,:)<=l))
            if ~good(u,n,a,l)
                good_a=false;
            end
        end
    end
end

n_config=size(data,1);



%%
% * Removing depth
idx=pos(:,2)<max_ze;
if any(idx==0)
    disp(['We removed data below the max depth ', num2str(max_ze), 'm which correspond to ', num2str(n_config-sum(idx)),' point(s)' ])
end
data=data(idx,:); pos=pos(idx,:);
n_config=size(data,1);

%%
% * Removing repetition
[~,idx,~]=unique(pos,'rows','first');
if length(idx)~=n_config
    disp(['Some configuration correspond to the same position, we keep the first configuration (which should minimized a and n and therefore the best). Removed data: ', num2str(length(idx)-n_config)])
end
data=data(idx,:); pos=pos(idx,:); k=k(idx);
n_config=size(data,1);


%%
% * Random sampling
if n_max<n_config
    disp(['Random sampling ',num2str(n_max),' among ',num2str(n_config),' points'])
    idx=datasample(1:n_config,n_max);
data=data(idx,:); pos=pos(idx,:);
n_config=size(data,1);
end



%%
% * PLOT
if plotit
    figure; plot(data,1:n_config,'x')
    set(gca, 'YDir', 'reverse'); set(gca,'xtick',[0:l])
    grid on; xlabel('Electrode position'); ylabel('Configuration')
    legend('A','B','M','N');
    
    figure;hold on
    plot(1:l,zeros(size(1:l)),'x'); plot(pos(:,1),pos(:,2),'o')
    xlabel('Electrode position');ylabel('depth')
    set(gca, 'YDir', 'reverse');
end


disp(['The number of configuration is ', num2str(n_config)])  
if n_config>4000
   pause()
end

end