function [ cam ] = getCam(varargin)
%GETCAM Function to establish figure and camera properties. Supporting
%function for plotMode and animateMode.
%
%
% Use like this:  (1) cam=getCam('load') --> define mat-file to load
%                 (2) cam=getCam('save') --> save current figure properties
%                     to mat file, define mat-file to save to
%                 (3) cam=getCam --> no save or load, simply get
%                     figure/camera properties from figure and store in cam
%                 (4) cam=getCam('C:\MyPath\To\Cam\File\cam.mat') loads cam
%                     from given file
%
%
% Knut Andreas Kvåle, 2014
%

if nargin==0
    disp('Adjust camera to set common settings for all plotting! Press any key to continue...')
    pause
    cam=getitall;
elseif nargin==1
    if ischar(varargin{1})
        if strcmp(varargin{1},'load')
            [filename,pathname]=uigetfile('*.mat','Load mat-file with figure properties...');
            cam=load([pathname,filename]);
            cam=cam.cam;
        elseif strcmp(varargin{1},'save')
            disp('Adjust camera to set common settings for all plotting! Press any key to continue...')
            pause
            cam=getitall;
            [filename,pathname]=uiputfile('*.mat','Save mat-file with figure properties...');
            save([pathname,filename],'cam');
        else
            cam=load(varargin{1});
            cam=cam.cam;
        end
    end
end

    function gotten=getitall()
    gotten.CameraPosition=get(gca,'CameraPosition');
    gotten.CameraTarget=get(gca,'CameraTarget');
    gotten.CameraUpVector=get(gca,'CameraUpVector');
    gotten.CameraViewAngle=get(gca,'CameraViewAngle');
    gotten.Position=get(gca,'Position');
    gotten.Projection=get(gca,'Projection');
    gotten.FigurePosition=get(gcf,'Position');
    gotten.FigurePosition=get(gcf,'Position');
    texthandle=findall(gca,'Type','text');
    gotten.TitlePosition=get(texthandle(1),'Position');
%     gotten.Light=get(findall(gcf,'type','light'));
    end
end

