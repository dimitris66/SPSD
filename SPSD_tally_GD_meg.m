% NOTE for Gio & Dimitrios, Mar 10, 2017, nima
% check this and set these valubles...
% maximafilt = 5 % number of passes of short moving avg for getting rid of small peaks
% selectmaximafilt = 10 % The number of passes of short moving avg for getting rid of small "selected" picks
% prctl_cut = 25 % The lower 25% of selected peaks will be discarded
% stdthresh = 1 % the lower percentage of the peaks that will be discarde. default=1 meaning that peaks smaller than stdev will be discarded
%--------------------section 1: find the peaks (local maxima)

%% Input:
% spsd_data: [timepoints]x[channels]
% type : string  e.g., 'EEG','MEG'
function [spindles] = SPSD_tally_GD_meg(spsd_data,type)

maximafilt = 5; % number of passes of short moving avg for getting rid of small peaks
selectmaximafilt = 10; % The number of passes of short moving avg for getting rid of small "selected" picks
prctl_cut = 25; % The lower 25% of selected peaks will be discarded
stdthresh = 2; % the lower percentage of the peaks that will be discarded. default=1 meaning that peaks smaller than stdev will be discarded

modalityname = type;
tsize = size(spsd_data,1);

%------------
display ('finding peaks...')
[aa bb] = size(spsd_data); % [timepoints]x[channels]
spsdmax.val=[];
spsdmax.ind=[];

for i =1:bb;
    [lmvall,indd] = lmax(spsd_data(:,i),maximafilt);
    spsdmax.val = [spsdmax.val lmvall];                            
    spsdmax.ind = [spsdmax.ind indd];                            
end
display ('finding peaks bigger than the threshold...')
spsdmax.z= find(spsdmax.val>(std2(spsd_data)*stdthresh)); % ****
% spsdmax.z= find(spsdmax.val>(std(spsdmax.val)*stdthresh)); % ****

% MEG
spsdmax.val=cell(1,bb);
spsdmax.ind=cell(1,bb);
for i =1:bb;
    [lmvall,indd] = lmax(spsd_data(:,i),maximafilt);
    spsdmax.val{i} = lmvall;                            
    spsdmax.ind{i} = indd;                            
end
display ('finding peaks bigger than the threshold...')

spsd_data_std = std(spsd_data);
spsdmax.z = cell(1,bb);
for i=1:bb
    spsdmax.z{i} = find(spsdmax.val{i}>(spsd_data_std(i)*stdthresh)); % indices of ("valid") maxima for every channel
end

%------------------section 2: 
% % % display('finding the number of peaks per sample point...')
% % % ipeak=[];
% % % for i=1:aa
% % %     tmp = [];
% % %     for j=1:bb
% % %          tmp = [tmp find(spsdmax.ind{j}(spsdmax.z{j}) == i)];
% % %     end
% % %     ipeak.ind{i} = tmp;
% % %     ipeak.size(i) = length(ipeak.ind{i});
% % % end    


display('finding the number of peaks per sample point...')
ipeak=[];
for i=1:aa
    tmp = [];
    for j=1:bb
        if ~isempty(find(spsdmax.ind{j}(spsdmax.z{j}) == i, 1))
            tmp = [tmp j];
        end
    end
    ipeak.ind{i} = tmp;
    ipeak.size(i) = length(ipeak.ind{i});
end    


%----------------plot
colordef none
figure;

%---plot spsd
display('plotting spsd...')
subplot(3,1,1)
plot(spsd_data,'y')
hold on
plot(mean(spsd_data,2),'r')
axis tight
xlabel('time')
ylabel('SPSD')
title(modalityname)

%-----plot spsd peak
display('plotting spsd peaks...')
subplot(3,1,2)
hold on
for i =1:bb;
    [lmvall,indd] = lmax(spsd_data(:,i),maximafilt);
    plot(indd(lmvall>(spsd_data_std(bb)*stdthresh)),lmvall(lmvall>(spsd_data_std(bb)*stdthresh)),'r*')  
%     plot(indd(lmvall>(std(lmvall)*stdthresh)),lmvall(lmvall>(std(lmvall)*stdthresh)),'r*')
end
%axis tight
xlabel('time')
title('SPSD peaks')

%------select spsd peaks based on number of active channels and plot it
display('selecting spsd peaks based on number of active channels...')
subplot(3,1,3)
hold on


tsize = ipeak.size.^2;
tsize(tsize < (std2(spsd_data)*stdthresh)) = 0;

% normalizing tsize
norm_tsize = tsize/max(tsize);
plot(norm_tsize)
axis tight
% finding the local maxima
[lmvall,indd] = lmax(norm_tsize,selectmaximafilt);
% trashing the lower percentage of the peaks (based on prctl_cut)
selectpeak=indd(lmvall>(prctile(lmvall,prctl_cut)));
ln = zeros(size(norm_tsize)); ln(indd(lmvall>(prctile(lmvall,prctl_cut)))) = lmvall(lmvall>(prctile(lmvall,prctl_cut)));
stem(ln,'LineStyle','-.','Color',[0 1 0],'MarkerEdgeColor','none','MarkerFaceColor','none')
plot(selectpeak,norm_tsize(selectpeak),'*g')
line([1 length(norm_tsize)], [prctile(lmvall,prctl_cut) prctile(lmvall,prctl_cut)],'Color',[1 0 0])
set(gca,'YTick',[])
xlabel('time')
title('spindle selection based on SPSD')
set(gcf,'InvertHardCopy','off')
%print('-dpdf',figname)
%close

%---------------
% display('Calculating number of sensors involved in significant peaks...')
% peakmtx =[];
% for i = 1:length(selectpeak)                 
%     peakmtx = [peakmtx;ipeak.size(:,selectpeak(i)-flank_win:selectpeak(i)+flank_win)];
% end
% peaksize = mean(peakmtx,2);
% 
% spindles = peakmtx;
spindles = [];


