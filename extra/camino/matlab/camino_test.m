addpath('../../../dtitools/matlab');
%% Test dt_track: a helicoid
if strcmp(questdlg('Run Tracking Tests?','Yes'),'Yes')
n1 = 64; n2 = 64; n3 = 64; N = 1000;
t = linspace(0,2*pi,N)';
c   = [n1/4*cos(t(:)) + n1/2,n2/4*sin(t(:)) + n2/2,n3/(2*pi)*t(:)] + 1;
c_t = [-n1/4*sin(t), n2/4*cos(t), n3/(2*pi)*ones(size(t))];
rc  = round(c); x = rc(:,1); y = rc(:,2); z = rc(:,3);
rc(x<1 | x>n1 | y<1 | y>n2 | z<1 | z>n3,:) = [];
D = zeros(n1,n2,n3,3,3); PD = zeros(n1,n2,n3,3);
for k = 1:length(rc)-1
  D(x(k),y(k),z(k),:,:) = c_t(k,:)'*c_t(k,:);
  PD(x(k),y(k),z(k),:)  = c_t(k,:)';
end;
FA=invariants(D,'fa'); ROI = double(FA>0.5);
figure;
disp('DT Tracking');
tic,tr_dt = dt_track(D,ROI); toc
streamline(tr_dt); view(3); axis equal; axis off
figure
disp('DT Tracking with noise');
D = D + 0.001*randn(size(D));
for x1 = 1:n1; for x2 = 1:n2; for x3 = 1:n3;
      D(x1,x2,x3,:,:) = sqrtm(squeeze(D(x1,x2,x3,:,:))'*squeeze(D(x1,x2,x3,:,:)));
end; end; end;
tic,tr_dt = dt_track(D,ROI); toc
streamline(tr_dt); view(3); axis equal; axis off
figure;
figure;disp('PD Tracking');
tic,tr_pd = pd_track(PD,ROI); toc
streamline(tr_pd); view(3); axis equal; axis off

%% Test datasynth
elseif strcmp(questdlg('Data Simulation','Yes'),'Yes')
disp('Data simulation');
V = datasynth(1,'../test/bmx7.scheme','brownian','N_walkers',10,'SNR',16,'seed',100000,'tauScale',1000,'qScale',1000);
end
%% REPORT ON TESTS: