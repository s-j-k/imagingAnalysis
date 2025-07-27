load('O:\sjk\DATA\imagingData\2p-opto\sk272\272-025-026\suite2p\plane0\Fall.mat')


Fnorm=F;
for ii = 1:size(F,1)
    Fnorm(ii,:)=(F(ii,:)-median(F(ii,:)))./median(F(ii,:));
end

figure; 
subplot(3,1,1);
imagesc(Fnorm(:,1:1000));hold on;
xline([24:100:1024],'w'); title('1000 frames after 10Hz');
clim([-0.15 1]);
subplot(3,1,2);imagesc(Fnorm(:,1001:2000));hold on;
xline([24:100:1024],'w'); title('1000-2000 frames after 10Hz');
clim([-0.15 1]);
subplot(3,1,3);imagesc(Fnorm(:,2001:3000));hold on;
xline([24:100:1024],'w'); title('2000-3000 frames after 10Hz');
clim([-0.15 1]);

figure;plot(mean(Fnorm,1));
ylim([-0.1 0.2]);
xline([20:100:1000],'r');title('1000 frames of white noise');


%%

load('O:\sjk\DATA\imagingData\2p-opto\sk272\272-040-041\suite2p\plane0\Fall.mat')

Fnorm=F;

for ii = 1:size(F,1)
    Fnorm(ii,:)=(F(ii,:)-median(F(ii,:)))./median(F(ii,:));
end

noLight=Fnorm(:,3084:4083);

figure; 
subplot(3,1,1);
imagesc(noLight(:,1:1000));hold on;
xline([24:100:1024],'w'); title('1000 frames of white noise');
clim([-0.15 1]);
subplot(3,1,2);imagesc(Fnorm(:,1:1000));hold on;
xline([24:100:1024],'w'); title('1-1000 frames after 10Hz');
clim([-0.15 1]);
subplot(3,1,3);imagesc(Fnorm(:,2001:3000));hold on;
xline([24:100:1024],'w'); title('1000-2000 frames after 10Hz');
clim([-0.15 1]);

figure;plot(mean(noLight,1));
ylim([-0.1 0.2]);
xline([20:100:1000],'r');title('1000 frames of white noise');