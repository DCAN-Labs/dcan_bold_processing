function varargout = motion_survival_plot(varargin)
% MOTION_SURVIVAL_PLOT MATLAB code for motion_survival_plot.fig
%      MOTION_SURVIVAL_PLOT, by itself, creates a new MOTION_SURVIVAL_PLOT or raises the existing
%      singleton*.
%
%      H = MOTION_SURVIVAL_PLOT returns the handle to a new MOTION_SURVIVAL_PLOT or the handle to
%      the existing singleton*.
%
%      MOTION_SURVIVAL_PLOT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MOTION_SURVIVAL_PLOT.M with the given input arguments.
%
%      MOTION_SURVIVAL_PLOT('Property','Value',...) creates a new MOTION_SURVIVAL_PLOT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before motion_survival_plot_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to motion_survival_plot_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help motion_survival_plot

% Last Modified by GUIDE v2.5 25-Jun-2014 15:27:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @motion_survival_plot_OpeningFcn, ...
                   'gui_OutputFcn',  @motion_survival_plot_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before motion_survival_plot is made visible.
function motion_survival_plot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to motion_survival_plot (see VARARGIN)

% Choose default command line output for motion_survival_plot
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes motion_survival_plot wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = motion_survival_plot_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function group_edit_Callback(hObject, eventdata, handles)
% hObject    handle to group_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of group_edit as text
%        str2double(get(hObject,'String')) returns contents of group_edit as a double



% --- Executes during object creation, after setting all properties.
function group_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to group_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function subject_edit_Callback(hObject, eventdata, handles)
% hObject    handle to subject_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of subject_edit as text
%        str2double(get(hObject,'String')) returns contents of subject_edit as a double


% --- Executes during object creation, after setting all properties.
function subject_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subject_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in group_browse_button.
function group_browse_button_Callback(hObject, eventdata, handles)
% hObject    handle to group_browse_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

group_folder_path = uigetdir(get(handles.group_edit,'String'),'Select group folder path');
set(handles.group_edit,'String',group_folder_path)

% --- Executes on button press in subject_browse_button.
function subject_browse_button_Callback(hObject, eventdata, handles)
% hObject    handle to subject_browse_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[subjlistfile,subjlistpath] = uigetfile([get(handles.group_edit,'String') filesep '*.txt'],'Select subject list file');
subject_list_file = [subjlistpath subjlistfile];
set(handles.subject_edit,'String',subject_list_file)


function fd_thresh_edit_Callback(hObject, eventdata, handles)
% hObject    handle to fd_thresh_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fd_thresh_edit as text
%        str2double(get(hObject,'String')) returns contents of fd_thresh_edit as a double

temp_value = str2double(get(hObject,'string'));
if isnan(temp_value)
      temp_value = 0.2; % set FD to the default
      msgbox(['The value you entered was not numeric.The FD threshold has been set to the default of ' num2str(temp_value) '.' 'Please enter a valid number if you would like to use a different FD threshold.'],'Invalid Entry')
end

set(hObject,'String',num2str(temp_value,'%1.2f'));

% --- Executes during object creation, after setting all properties.
function fd_thresh_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fd_thresh_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dvar_thresh_edit_Callback(hObject, eventdata, handles)
% hObject    handle to dvar_thresh_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dvar_thresh_edit as text
%        str2double(get(hObject,'String')) returns contents of dvar_thresh_edit as a double

temp_value = round(str2double(get(hObject,'string')));
if isnan(temp_value)
      temp_value = 20; % set DVAR to the default
      msgbox(['The value you entered was not numeric.The DVAR threshold has been set to the default of ' num2str(temp_value) '.' 'Please enter a valid number if you would like to use a different DVAR threshold.'],'Invalid Entry')
end

set(hObject,'String',num2str(temp_value,'%d'));

% --- Executes during object creation, after setting all properties.
function dvar_thresh_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dvar_thresh_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot_button.
function plot_button_Callback(hObject, eventdata, handles)
% hObject    handle to plot_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

group_folder_path = get(handles.group_edit,'String');
subject_list_file = get(handles.subject_edit,'String');
FD_threshold = str2double(get(handles.fd_thresh_edit,'String'));
DVAR_threshold = str2double(get(handles.dvar_thresh_edit,'String'));

[subject_list_cell,subject_survival] = power_2014_motion_survival_plot(group_folder_path,subject_list_file,FD_threshold);

imagesc(subject_survival)
set(handles.survival_plot,'xtick',[],'ytick',1:length(subject_list_cell),'yticklabel',subject_list_cell);
ylabel('Subject Visits')
title('Remaining time in minutes after FD thresholding')
colorbar
colormap('hot')
