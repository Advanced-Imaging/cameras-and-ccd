%  Tempering project code made by Nir Goldfriend of Ori Katz Lab 

clear 
close all


inputByCam = 1; % 0 for given video, 1 for live video
ifStop = 0; % stop after the end of the video
ifshow = 1; % optional not to show the figures
ifshowHIST=0;
if inputByCam == 1
% end
try
    delete(cam_direct)
    delete(cam_split)
    clear cam_split cam_direct
catch
    
    imaqreset
end
end
numFrames = 100; %length of the video
axes_bw=1;
if axes_bw ~=1 && axes_bw~= 0  
    axes_bw =0;
end

CrossCorr = zeros(2,numFrames); % maximum CC (CrossCorrelation) normalized
CrossCorrLast = zeros(2,numFrames);
CrossCorr(:,1) = [1 ;1] ;
CrossCorrLast(:,1)  = [1 ;1] ;
x = zeros(1,numFrames); % creating arrays for spatial representations 
y = zeros(1,numFrames);  % creating arrays for spatial representations
z = zeros(1,numFrames);
t = zeros(1,numFrames);  % creating arrays for time

 figure_split=figure;
 figure_direct=figure;
 figure_live=figure;
 figure_split.NumberTitle='off';
 figure_direct.NumberTitle='off';
 figure_direct.Position=[22 576 560 415];
 figure_split.Position=[44 98 560 420];
 figure_live.NumberTitle='off';
 figure_live.Color=axes_bw.*[1 1 1];
 figure_live.Name='Live View';
 if ifshowHIST==1
    figure_hist=figure;
 end
 ExposureTime=15000; %ExposureTime for the cameras
 Timeout=50; %Timeout for the cameras
CrossCorr(:,1) = [1;1] ; % The CrossCorr for self image t=0
CrossCorrLast(:,1)  = [1;1] ; % The CrossCorrLast for for self image t=0
%  Motorz = initMotor(27000356); % The motor for z axis
%  Motorx = initMotor(27000357); % The motor for x axis

%%
counter = 1;
if inputByCam==1
    imaqreset
    
   
try
    cam_direct = videoinput('gige', 1, 'Mono8');
    cam_split = videoinput('gentl', 1, 'Mono8');  
catch
        
    cam_split = videoinput('gentl', 2, 'Mono8');
    cam_direct = videoinput('gige', 1, 'Mono8');

     
end
%     triggerconfig(cam_split,'manual');
%     triggerconfig(cam_direct,'manual');
    cam_direct.Timeout=Timeout;
    cam_split.Timeout=Timeout;
    
    src_split = getselectedsource(cam_split);
    src_direct = getselectedsource(cam_direct);
    src_split.Gain=0;
    src_direct.Gain=0;

     exposure_split=(17/15)*ExposureTime;
    exposure_direct=(20/15)*ExposureTime;
     src_split.ExposureTime=exposure_split;
     src_direct.ExposureTimeAbs=exposure_direct;

    
    
%     start(cam_split)
%     start(cam_direct)
%         trigger(cam_split);trigger(cam_direct)
      pic_split = getsnapshot(cam_split);
      pic_direct= getsnapshot(cam_direct);
     

    pic_1_split = double(pic_split(:,:,:,1)); %first pic
    pic_1_direct = double(pic_direct(:,:,:,1)); %first pic
%     
%     switch mean(max(pic_1_split))
%         case  mean(max(pic_1_split))< 100
%         src_split.ExposureTime=exposure_split*30;
%         pic_split = getsnapshot(cam_split);
%         pic_1_split = double(pic_split(:,:,:,1)); %retake first pic
%         case  mean(max(pic_1_split))> 120 && mean(max(pic_1_split)) <252
%         src_split.ExposureTime=exposure_split/12;
%         pic_split = getsnapshot(cam_split);
%         pic_1_split = double(pic_split(:,:,:,1)); %retake first pic
%         case  mean(max(pic_1_split))> 120 && mean(max(pic_1_split)) >252
%         src_split.ExposureTime=exposure_split/120;
%         pic_split = getsnapshot(cam_split);
%         pic_1_split = double(pic_split(:,:,:,1)); %retake first pic
%       
%     end
%     
%     
%     switch mean(max(pic_1_direct))
%         case  mean(max(pic_1_direct))< 100
%         src_direct.ExposureTime=exposure_direct*30;
%         pic_direct = getsnapshot(cam_direct);
%         pic_1_direct = double(pic_direct(:,:,:,1)); %retake first pic
%         case  mean(max(pic_1_direct))> 120 && mean(max(pic_1_direct)) <252
%         src_direct.ExposureTime=exposure_direct/12;
%         pic_direct = getsnapshot(cam_direct);
%         pic_1_direct = double(pic_direct(:,:,:,1)); %retake first pic
%         case  mean(max(pic_1_direct))> 120 && mean(max(pic_1_direct)) >252
%         src_direct.ExposureTime=exposure_direct/120;
%         pic_direct = getsnapshot(cam_direct);
%         pic_1_direct = double(pic_direct(:,:,:,1)); %retake first pic
%       
%     end
%     
%         

     
    
    lastPic_split=pic_1_split;
    lastPic_direct=pic_1_direct;
    vid_split = zeros(size(pic_1_split,1),size(pic_1_split,2),numFrames); %keep a new video
    vid_direct = zeros(size(pic_1_direct,1),size(pic_1_direct,2),numFrames); %keep a new video
end

%%
if inputByCam==1
    tic


    while counter<=numFrames %true%
        

     if counter>20
         pause (2800)
     end
    try
        pic_split = getsnapshot(cam_split);
        pic_direct = getsnapshot(cam_direct);
        pic_split=double(pic_split(:,:,:,1));
        pic_direct=double(pic_direct(:,:,:,1));
    catch

         imaqreset
               cam_split = videoinput('gentl', 2, 'Mono8');
                cam_direct = videoinput('gige', 1, 'Mono8');
                cam_direct.Timeout=Timeout+10;
                cam_split.Timeout=Timeout+10;

                src_split = getselectedsource(cam_split);
                src_direct = getselectedsource(cam_direct);

                src_split.Gain=0;
                src_direct.Gain=0;

                exposure_split=(12/15)*ExposureTime;
                exposure_direct=(14/15)*ExposureTime;
                src_split.ExposureTime=exposure_split;
                 src_direct.ExposureTimeAbs=exposure_direct;


               
        pic_split = getsnapshot(cam_split);
        pic_direct = getsnapshot(cam_direct);
        pic_split=double(pic_split(:,:,:,1));
        pic_direct=double(pic_direct(:,:,:,1));
        
    end
       

        if ifshowHIST==1
        figure(figure_hist)
        [intensity_split,counts_s ] = (histOfPic717(pic_split,0)); % histogram of intensity
        [intensity_direct,counts_d ] = (histOfPic717(pic_direct,0)); % histogram of intensity
%       semilogy(intensity,counts,'.r','MarkerSize',5);
        xlabel('Brightness[AU for 8bit(256) greyscale method]','FontSize',10)
        ylabel('Histogram Counts','FontSize',14)
        hold on
        grid on
        end
        
       t(counter) = toc;
        [CrossCorrTemp_split,xTemp_split,yTemp_split] = CrossCorrelationNorm717(pic_1_split,pic_split,1,figure_split); %find pic in pic1
        [CrossCorrTemp_direct,xTemp_direct,yTemp_direct] = CrossCorrelationNorm717(pic_1_direct,pic_direct,1,figure_direct); %find pic in pic1
          
            CrossCorrTempLast_split = CrossCorrelationNorm717(lastPic_split,pic_split,0);
            CrossCorrTempLast_direct = CrossCorrelationNorm717(lastPic_direct,pic_direct,0);
            CrossCorrLast(:,counter) = [CrossCorrTempLast_split;CrossCorrTempLast_direct]; % 1st_row  is the split and the 2nd is the direct please follow
            lastPic_split = pic_split;
            lastPic_direct = pic_direct;

        
        F_direct = fftshift(fft2(pic_direct-mean(pic_direct(:))));  
        F_split = fftshift(fft2(pic_split-mean(pic_split(:)))); %it can be helpful to shift the zero-frequency components to the center
        CrossCorr(:,counter) = [CrossCorrTemp_split;CrossCorrTemp_direct];
    
            x(counter) = xTemp_split;
            y(counter) = yTemp_split;
            vid_split(:,:,counter) = pic_split;
            vid_direct(:,:,counter) = pic_direct;
%                 if exist('TempFile_.mat', 'file')==2
%                     delete('TempFile_.mat');
%                 end
             save(['TempFile_',int2str(counter)],'CrossCorr','CrossCorrLast','pic_split','pic_direct','t','x','y')
     pause(0.1);

                       
         if ifshow==1
            figure(figure_live);
            set(gcf,'position', [852 49 1066 948]);
            %figure_live.CurrentAxes.Color='w';
            colormap hot
            subplot(3,3,2)
            imagesc(pic_split(1:300,1:300)) % show the picture
            colorbar
            title(['Picture of Split Beam ',num2str(counter)],'FontSize',16,'Color', (1-axes_bw).*[1 1 1])
            axis image
            ax = gca;
            ax.XColor=(1-axes_bw).*[1 1 1]; ax.YColor=(1-axes_bw).*[1 1 1];
            
            
            subplot(3,3,3)
            imagesc(pic_direct(1:300,1:300)) % show the picture  
            colorbar
            axis image
             ax = gca;
            ax.XColor=(1-axes_bw).*[1 1 1]; ax.YColor=(1-axes_bw).*[1 1 1];
            drawnow
            title(['Picture of Direct Beam ',num2str(counter)],'FontSize',16,'Color', (1-axes_bw).*[1 1 1])
            subplot(3,3,1)
             d_x=x(1:counter-1); %
            d_y=y(1:counter-1); %
            d_r=(sqrt(d_x.^2+d_y.^2)); %
            
            plot(t(1:counter-1),[d_r ; d_x ; d_y]);
            legend('\Delta r ','\Delta x ','\Delta y ')
            clear d_r d_x d_y
            xlabel('time [s]','Color', (1-axes_bw).*[1 1 1]);
            ylabel('distance of max [frames]','Color', (1-axes_bw).*[1 1 1]);
            ax = gca;
            ax.XColor=(1-axes_bw).*[1 1 1]; ax.YColor=(1-axes_bw).*[1 1 1];
            title('Translation','FontSize',12,'Color', (1-axes_bw).*[1 1 1]);
            subplot(3,3,4)
            imagesc(log10((abs(F_split))),max(log10((abs(F_split(:)))))+[-3 0]);
            colorbar
            title('Fourier of Intensity Splited Beam, log10 Scale','Color', (1-axes_bw).*[1 1 1])
            axis image
            ax = gca;
            ax.XColor=(1-axes_bw).*[1 1 1]; ax.YColor=(1-axes_bw).*[1 1 1];
            subplot(3,3,7)
            imagesc(log10((abs(F_direct))),max(log10((abs(F_direct(:)))))+[-3 0]);
            colorbar
            title('Fourier of Intensity Direct Beam, log10 Scale','Color', (1-axes_bw).*[1 1 1])
             axis image
            ax = gca;
            ax.XColor=(1-axes_bw).*[1 1 1]; ax.YColor=(1-axes_bw).*[1 1 1];
            
            subplot(3,3,[5 6 8 9]);
            plot(t(1:counter-1),CrossCorr(:,1:counter-1),'Linewidth',2)
            xlabel('Time [s]','Color', (1-axes_bw).*[1 1 1]);
            ylabel('CrossCorrelation for both beams','Color', (1-axes_bw).*[1 1 1]);
            legend('Split Beam','Direct Beam')
 
            ax.XColor=(1-axes_bw).*[1 1 1]; ax.YColor=(1-axes_bw).*[1 1 1];           
            
         end
     if counter>10 && counter<31 
    moveStep(0.0005,Motorx)
    x(counter)=0.5*(counter-10); % in mm
    end
    x(31:32)=x(30).*ones(1,2);
    if counter>32 && counter<52
            moveStep(-0.0005,Motorx)
    x(counter)=x(32)-0.5*(counter-32); % in mm
    end
    
     if counter>56 && counter<70 
    moveStep(0.0005,Motorz)
    z(counter)=0.5*(counter-56); % in mm
    end
    z(31:32)=z(30).*ones(1,2);
    if counter>73 && counter<100
            moveStep(-0.0005,Motorz)
    z(counter)=z(73)-0.5*(counter-73); % in mm
    end
   
      counter = counter+1;   
    end
% stop(cam_direct)
% stop(cam_split)
delete(cam_direct)
delete(cam_split)
clear cam_split cam_direct
%%
else
    pic_1_split = double(vid_split(:,:,1)); %first pic
    pic_1_direct = double(vid_direct(:,:,1)); %first pic
    lastPic_split=pic_1_split;
    lastPic_direct=pic_1_direct;
    ifshow = 1;
    numFrames = min(numFrames,size(vid_direct,3)); %don't go over the size
    counter = 0;
    
    for counter=2:numFrames %true%
        pic_split=double(vid_split(:,:,counter)); %equivalent to g snapshot
        pic_direct=double(vid_direct(:,:,counter)); %equivalent to g snapshot
        
        %whos pic_1 pic
        % get(gcf,'position')
        if ifshowHIST==1
         figure(figure_hist)
        [intensity_split,counts_s ] = (histOfPic717(pic_split,0)); % histogram of intensity
        [intensity_direct,counts_d ] = (histOfPic717(pic_direct,0)); % histogram of intensity
%       semilogy(intensity,counts,'.r','MarkerSize',5);
        xlabel('Brightness[AU for 8bit(256) greyscale method]','FontSize',10)
        ylabel('Histogram Counts','FontSize',14)
        hold on
        grid on
        end
        
        [CrossCorrTemp_split,xTemp_split,yTemp_split] = CrossCorrelationNorm717(pic_1_split,pic_split,1,figure_split); %find pic in pic1
        [CrossCorrTemp_direct,xTemp_direct,yTemp_direct] = CrossCorrelationNorm717(pic_1_direct,pic_direct,1,figure_direct); %find pic in pic1
        if counter>2
            
            CrossCorrTempLast_split = CrossCorrelationNorm717(lastPic_split,pic_split,0);
            CrossCorrTempLast_direct = CrossCorrelationNorm717(lastPic_split,pic_split,0);
            CrossCorrLast(:,counter) = [CrossCorrTempLast_split;CrossCorrTempLast_direct]; % 1st_row  is the split and the 2nd is the direct please follow
            lastPic_split = pic_split;
            lastPic_direct = pic_direct;
        end
        
        F_direct = fftshift(fft2(pic_direct-mean(pic_direct(:))));  
        F_split = fftshift(fft2(pic_split-mean(pic_split(:)))); %it can be helpful to shift the zero-frequency components to the center
        CrossCorr(:,counter) = [CrossCorrTemp_split;CrossCorrTemp_direct];
    
            x(counter) = xTemp_split;
            y(counter) = yTemp_split;
            
%                 if exist('TempFile_.mat', 'file')==2
%                     delete('TempFile_.mat');
%                 end
                                   
         if ifshow==1
            figure(figure_live);
            set(gcf,'position', [852 49 1066 948]);
            %figure_live.CurrentAxes.Color='w';
            colormap hot
            subplot(3,3,2)
            imagesc(pic_split(1:300,1:300)) % show the picture
            colorbar
            title(['Picture of Split Beam ',num2str(counter)],'FontSize',16,'Color', 'k')
            axis image
%             ax = gca;
%             ax.XColor='w'; ax.YColor='w';
            subplot(3,3,3)
            imagesc(pic_direct(1:300,1:300)) % show the picture
            colorbar
            axis image
%             ax = gca;
%             ax.XColor='w'; ax.YColor='w';
            drawnow
            title(['Picture of Direct Beam ',num2str(counter)],'FontSize',16,'Color', 'k')
            subplot(3,3,1)
             d_x=x(1:counter-1); %
            d_y=y(1:counter-1); %
            d_r=(sqrt(d_x.^2+d_y.^2)); %
            plot(t(1:counter-1),[d_r ; d_x ; d_y]);
            legend('\Delta r ','\Delta x ','\Delta y ')
            clear d_r d_x d_y
            xlabel('time [s]','Color', 'k');
            ylabel('distance of max [frames]','Color', 'k');
%             ax = gca;
%             ax.XColor='w'; ax.YColor='w';
            title('Translation','FontSize',12,'Color', 'k');
            subplot(3,3,4)
            imagesc(log10((abs(F_split))),max(log10((abs(F_split(:)))))+[-3 0]);
            colorbar
            title('Fourier of Intensity Splited Beam, log10 Scale','Color', 'k')
            axis image
%             ax = gca;
%             ax.XColor='w'; ax.YColor='w';
            subplot(3,3,7)
            imagesc(log10((abs(F_direct))),max(log10((abs(F_direct(:)))))+[-3 0]);
            colorbar
            title('Fourier of Intensity Direct Beam, log10 Scale','Color', 'k')
            axis image
%             ax = gca;
%             ax.XColor='w'; ax.YColor='w';
            subplot(3,3,[5 6 8 9]);
            plot(t(1:counter-1),CrossCorr(:,1:counter-1),'Linewidth',2)
            xlabel('Time [s]','Color', 'k');
            ylabel('CrossCorrelation for both beams','Color', 'k');
            legend('Split Beam','Direct Beam')
%             ax = gca;
%             ax.XColor='w'; ax.YColor='w';
 
            
         end
    end
    
%              implay(givenVid);


end
%%

hold off
 if counter>=30 && inputByCam==1
     
     save(datestr(now,30),'CrossCorrTempLast_direct','CrossCorrTempLast_split','CrossCorrTemp_direct'...
     ,'CrossCorrTemp_split','vid_direct',...                                
    'vid_split','t','x','y');
               
     savefig(figure_direct,['figure_direct',datestr(now,30)],'compact');    
      savefig(figure_hist,['figure_hist',datestr(now,30)],'compact');    
       savefig(figure_live,['figure_live',datestr(now,30)],'compact');    
        savefig(figure_split,['figure_split',datestr(now,30)],'compact');    
      for n=1:99
     delete (['TempFile_',int2str(n),'.mat'])
     end
 end