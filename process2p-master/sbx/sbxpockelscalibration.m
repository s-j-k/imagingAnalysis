
% Pockels calibration

global scanbox_h sbconfig

pockels_val = round(linspace(0,255,32));
pockels_cal = zeros(size(pockels_val));

sb_pockels_lut_identity;    % load identity
sb_setparam(512,0,4);
sb_scan;

fig = figure('MenuBar','None','NumberTitle','off','Name','Pockels Calibration');
xlabel('DAC Value');
ylabel('Power (W)');
xlim([0 256]);
hold on;
for i = 1:length(pockels_val)
    sb_pockels(0,pockels_val(i));
    pause(5);
    pockels_cal(i) = powermeter_read;
    plot(pockels_val(i),pockels_cal(i),'b.','markersize',20);
end

sb_abort

% find max and fit curve

[~,idx] = max(pockels_cal);

xdata = pockels_val(1:idx);
ydata = pockels_cal(1:idx);
p = lsqcurvefit(@(x,xdata) x(1) - x(2) * cos(xdata/x(3)*2*pi),[mean(pockels_cal) mean(pockels_cal)*0.5 8*idx],xdata,ydata);
plot(xdata,p(1)-p(2) * cos(xdata*2*pi/p(3)),'r-')

yy = linspace(p(1)-p(2),p(1)-p(2)*cos(xdata(idx)*2*pi/p(3)),256);
pockels_lut = uint8(round(acos(-(yy-p(1))/p(2))*p(3)/pi/2));

% save 

sbxroot = fileparts(which('scanbox'));
save([sbxroot '\pockelscal.mat'],'pockels_lut'); % save otcal in the core directory (otcalibration is too big!)

% validation
sbconfig.pockels_lut = pockels_lut;
for(i=1:256)
    sb_pockels_lut(i,sbconfig.pockels_lut(i));
end

sb_setparam(512,0,4);
sb_scan;

fig = figure('MenuBar','None','NumberTitle','off','Name','Pockels Calibration Validation');
xlabel('DAC Value');
ylabel('Power (W)');
xlim([0 256]);
hold on;
for i = 1:length(pockels_val)
    sb_pockels(0,pockels_val(i));
    pause(5);
    pockels_cal(i) = powermeter_read;
    plot(pockels_val(i),pockels_cal(i),'b.','markersize',20);
end

sb_abort




% 
% 
% 
% sbconfig.optocal = [];          % make it empty
% 
% h = guihandles(scanbox_h);
% 
% oval = sbconfig.optoval;        % sequence of current values
% 
% % set animal and unit numbers
% 
% h.animal.String = 'xx0'; h.animal.Callback(h.animal,[]);
% h.unit.String = '000';   h.unit.Callback(h.unit,[]);
% h.expt.String = '000';   h.expt.Callback(h.unit,[]);
% 
% cd([h.dirname.String '\' h.animal.String]);
% 
% system('del xx0*');                                 % removes xx0 files
% 
% tri_send('KBY',0,12,0); % set knobby to super fine
% 
% h.tfilter.Value = 3; h.tfilter.Callback(h.tfilter,[]);  % accumulate
% 
% % set knobby params + table
% 
% h.zrange.String = num2str(sbconfig.optorange);  
% h.zstep.String = num2str(sbconfig.optostep);
% h.framesperstep.String = num2str(sbconfig.optoframes);
% h.zrange.Callback(h.zrange,h);   % update fields and table
% 
% % set number of frames to collect
% 
% h.frames.String = num2str(h.knobby_table.Data(end,end) + str2double(h.framesperstep.String));
% h.frames.Callback(h.frames,[]);  % total number of frames to collect
% 
% % Get the button
% g = findobj(scanbox_h,'tag','grabb');
% 
% % assumes pollen grains in focus
% 
% tri_send('KBY',0,0,40); % move up by 35um
% pause(0.2);
% tri_send('KBY',0,30,0); % zero knobby
% pause(0.2);
% 
% % make sure it returns to first location
% 
% h.returnbox.Value=1; h.returnbox.Callback(h.returnbox,[]);
% 
% % collect the data
% for i = 1:length(oval)
%     h.knobby_enable.Value = 1; h.knobby_enable.Callback(h.knobby_enable,[]);
%     % set the optoslider
%     h.optoslider.Value = oval(i);
%     h.optoslider.Callback(h.optoslider,h);
%     h.grabb.Callback(h.grabb,[]);
%     pause(2);
% end
% 
% % process
% 
% clear otcalibration;
% 
% A = sbxread('xx0_000_000',0,1);
% A = squeeze(A(1,:,:,:));
% N = min(size(A))-80;    % leave out margin an focus on center
% yidx = round(size(A,1)/2)-N/2 + 1 : round(size(A,1)/2)+ N/2;
% xidx = round(size(A,2)/2)-N/2 + 1 : round(size(A,2)/2)+ N/2;
% 
% wb = waitbar(0,'Processing...  Please wait.');
% for i = 1:length(oval)
%     fn = sprintf('xx0_000_%03d',i-1);
%     z = sbxreadzstack(fn);    
%     m = squeeze(mean(z,3));
%     otcalibration(i).m = m;
%     [u,v] = fftalign(otcalibration(i).m,otcalibration(1).m);
%     otcalibration(i).uv = [u v];
%     z2 = circshift(z,[u v 0]); 
%     otcalibration(i).z = z2; 
%     xc = zeros(1,size(z,3));
%     for k = 0:size(z,3)-1
%         zs = circshift(otcalibration(i).z,[0 0 k]);
%         A = zs(yidx,xidx,:); B = otcalibration(1).z(yidx,xidx,:);
%         [rho,p] = corrcoef(A(:),B(:));
%         xc(k+1) = rho(1,2);
%     end
%     otcalibration(i).xc = xc;
%     waitbar(i/length(oval),wb);
% end
% delete(wb);
% global info
% fclose(info.fid);
% info = [];
% 
% % compute peak locations
% th = linspace(0,2*pi,size(z,3)+1);
% th = th(1:end-1);
% 
% for i = 1:length(oval)
%     otcalibration(i).th = angle(sum(exp(1i*th).*otcalibration(i).xc));
% end
% 
% depth = unwrap([otcalibration(:).th]) * abs(sbconfig.optorange) /(2*pi);
% 
% f = figure;
% plot(oval,-depth,'-o');
% otcoeff = polyfit(oval,-depth,2);
% hold on
% plot(0:1700,polyval(otcoeff,0:1700),'r-');
% xlabel('DAC ETL Value');
% ylabel('Depth (um)');
% title('Optotune Calibration');
% f.NumberTitle = 'off';
% f.MenuBar = 'none';
% 
% sbxroot = fileparts(which('scanbox'));
% save([sbxroot '\otcal.mat'],'otcoeff'); % save otcal in the core directory (otcalibration is too big!)
% clear otcalibration                     % delete huge variable
% 
% 
