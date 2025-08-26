%% Enhance reproducability!
%% Mirko Articus, October 2022
% This function is meant to document for all plots, where the shown data came
% from, which units were involved, which script was used etc...
% The original fullpath, scriptname, date and some more info will be printed in
% a textbox on the plot. All given information will additionally be saved in a struct.
% You should place this function right before the figure is saved/exported.

%% Input
% fig - handle to the figure you are about to annotate and save
% fullpath - where the figure is about to be saved
% scriptnName - dame of the script used to produce this figure
% to get scriptName use "mfilename('fullpath')" as input when calling the
% function
% saveFig - logical indicating if figure is saved and only if the fullpath
% will be printed and the info file will be created

% varargin - varname1,var1,varname2,var2...
% possible inputs ...
% dName - name of used d-struct
% uidx - indices for units included to produce this plot

% inputs printed on the plot
% always: date, fullpath, scriptName
% info about datasource -> names/fullpaths of mat-files: in varargin
% preceeded by 'dSrc1'/.../'dSrc4' -> 4 inputs allowed

%% Output
% annot - handle to the annotation textbox
% srcInfo - struct containing info plotted in textbox and some more. The
% file is saved as "info_<filename>.mat" on the same path as the plot
% itself

%%
function [annot, srcInfo] = docDataSrc(fig,fullpath,scriptName,saveFig,varargin)
%% textbox
annotStr = {date; fullpath;['Used Script: ' scriptName '.m']};
if ~saveFig; annotStr{2,1}='notsaved';end

% define additional info you'd like to be plotted on the figure
varNamePos = 1:2:numel(varargin);
addStr = varNamePos(contains(varargin(varNamePos),{'dSrc1','dSrc2','dSrc3','dSrc4',}))+1;
for asx =1:numel(addStr)
    annotStr = [annotStr;varargin(addStr(asx))];
end

spacerX =.35;
fig.Position = fig.Position + [0 spacerX*numel(annotStr) 0 spacerX*numel(annotStr)];
annotPos = [0.01 1 .1 .0]; % upleft corner


annot = annotation('textbox',annotPos,'String',annotStr,'FitBoxToText','on',...
    'Interpreter','none','EdgeColor','none','FontSize',6);
drawnow % helps if matlab is going "too fast

%% create and save info mat-file
if saveFig
    srcInfo.scriptName = scriptName;
    srcInfo.fullpath = fullpath;
    for vx = 1:2:numel(varargin)
        srcInfo.(varargin{vx}) = varargin{vx+1};
    end
    
    [filepath,name,~] = fileparts(fullpath);
    save([filepath filesep 'info_' name '.mat'],'srcInfo')
end
end