function [data,pos,n_config,k]=configuration(method, elec_n, depth_max, config_max, plotit)
%% CONFIGURATION compute the different configuration possible for a total
% number of electrode. There is no unit, you can think of it as an index
% for the matrix.
%
% INPUT:
% * method      defining the method of configuration 'wenner', 'pole-pole',
%               'dipole-dipole',...
% * elec_n      total number of electrode [no unit]
% * depth_max   max depth. Above this depth, all configuration will be
%               removed
% * config_max  Total number of configuration. If the code generate more
%               than config_max, configuration will be remove by how close it
%               is from its neighbourhood configuration up to config_max
%               configurations.
% * plotit      Bolean (1,0) decide for ploting or not
%
% OUPUT:
% * data        in coelec_numn order:
%               1: +ve emmiteur electrode
%               2: -ve emmiteur electrode
%               3: +ve receiveur electrode
%               4: -ve receiveur electrode
% * pos         average position of the measure for an homogenous media.
%               This is compute in x as the middle of the electrodes and z
%               as the investigation depth
% * n_config    finaelec_n total number of configuration
% * k           

assert(elec_n>4 && mod(elec_n,1)==0,'Invalid elec_n value')
assert(mod(depth_max,1)==0,'Invalid depth_max value')
assert(mod(config_max,1)==0,'Invalid config_max value')
assert(plotit==1 || plotit==0,'Invalid plotit value')

%%
% * Defining for the seelec_nected method
switch method
    case 'wenner'
        A_fx    = @(u,n,a) u;
        M_fx    = @(u,n,a) u+a;
        N_fx    = @(u,n,a) u+a*2;
        B_fx    = @(u,n,a) u+a*3;
        good    = @(u,n,a,elec_n) u+(a+1)*3<=elec_n && n<=1;
        ze_fx   = @(u,n,a) 0.173*a*3;
        xa_fx   = @(u,n,a) u+a*1.5;
        k_fx    = @(n,a) 2*pi*a;
    case 'poelec_ne-poelec_ne'
        A_fx    = @(u,n,a) u;
        M_fx    = @(u,n,a) u+a;
        N_fx    = @(u,n,a) NaN;
        B_fx    = @(u,n,a) NaN;
        good    = @(u,n,a,elec_n) u+(a+1)<=elec_n && n<=1;
        ze_fx   = @(u,n,a) 0.35*a;
        xa_fx   = @(u,n,a) u+a*0.5;
        k_fx    = @(n,a) 2*pi*a;
    case 'wenner-schelec_numberger'
        zec=[.173 .186 .189 .190];
        A_fx    = @(u,n,a) u;
        M_fx    = @(u,n,a) u+a*n;
        N_fx    = @(u,n,a) u+a*n+a;
        B_fx    = @(u,n,a) u+a*n+a+a*n;
        good    = @(u,n,a,elec_n) (u+(a+1)*n+(a+1)+(a+1)*n)<=elec_n;
        ze_fx   = @(u,n,a) zec(min(n,length(zec)))*(a*n+a+a*n);
        xa_fx   = @(u,n,a) u+n*a+a*0.5;
        k_fx    = @(n,a) pi*n*(n+1)*a;
    case 'dipole-dipole'
        zec=[.139 .174 .192 .203 .211 .216 .22 .224 .255];
        B_fx    = @(u,n,a) u;
        A_fx    = @(u,n,a) u+a;
        M_fx    = @(u,n,a) u+a+n*a;
        N_fx    = @(u,n,a) u+a+a*n+a;
        good    = @(u,n,a,elec_n) (u+(a+1)+(a+1)*n+(a+1))<=elec_n;
        ze_fx   = @(u,n,a) zec(min(n,length(zec)))*(a+a*n+a);
        xa_fx   = @(u,n,a) u+a+n*a*0.5;
        k_fx    = @(n,a)   pi*n*(n+1)*(n+2)*a;
    otherwise
        error('unvalid variabel method')
end


%%
% Computing all possibel configuration
config_n_pos = ceil(elec_n^3/3)/2;
data        = nan(config_n_pos,4);
pos         = nan(config_n_pos,2);
k           = nan(config_n_pos);
i=0;
u=0;
a=0;
good_u=true;
while good_u
    u=u+1;
    n=0;
    if ~good(u,n,a,elec_n) || ze_fx(u,n+1,a)> depth_max
        good_u=false;
    else
        good_n=true;
    end
    while good_n
        n=n+1;
        a=0;
        if ~good(u,n,a,elec_n) || ze_fx(u,n,a+1)> depth_max
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
            assert(any(data(i,:)<=elec_n))
            if ~good(u,n,a,elec_n) || ze_fx(u,n,a+1)> depth_max
                good_a=false;
            end
        end
    end
end
n_config=size(data,1);


%%
% * Removing depth
idx=pos(:,2)<depth_max;
if any(idx==0) && plotit
    disp(['We removed data beelec_now the max depth ', num2str(depth_max), 'm which correspond to ', num2str(n_config-sum(idx)),' point(s)' ])
end
data=data(idx,:); pos=pos(idx,:);
n_config=size(data,1);


%%
% * Remove by clustering...

if config_max<n_config
    idx = kmeans(pos,config_max);
    data_n = nan(config_max,4);
    pos_n  = nan(config_max,2);
    k_n    = nan(config_max);
    for i=1:config_max
        u=find(idx==i);
        [~, idx2]=max(k(u));
        data_n(i,:)=data(u(idx2),:); pos_n(i,:)=pos(u(idx2),:); k_n(i)=k(u(idx2));
    end
    data=data_n;
    pos=pos_n;
    k=k_n;
    n_config=config_max;
end


%%
% * PlOT
if plotit
    figure; plot(data,1:n_config,'x')
    set(gca, 'YDir', 'reverse'); set(gca,'xtick',0:elec_n)
    grid on; xlabel('electrode position'); ylabel('Configuration')
    legend('A','B','M','N');
    
    figure;hold on
    plot(1:elec_n,zeros(size(1:elec_n)),'x'); plot(pos(:,1),pos(:,2),'o')
    xlabel('electrode position');ylabel('depth')
    set(gca, 'YDir', 'reverse');
    
    disp(['The number of configuration is ', num2str(n_config)])
    if n_config>4000
        pause()
    end
end

end