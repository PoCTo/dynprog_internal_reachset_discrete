function varargout = gui(varargin)
% GUI MATLAB code for gui.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui

% Last Modified by GUIDE v2.5 16-Dec-2012 05:47:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_OutputFcn, ...
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

global opts;
opts = struct('ellipses','','approximations','','centers','','directions','');
% End initialization code - DO NOT EDIT


% --- Executes just before gui is made visible.
function gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui (see VARARGIN)

% Choose default command line output for gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in menuSolver.
function menuSolver_Callback(hObject, eventdata, handles)
% hObject    handle to menuSolver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns menuSolver contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menuSolver


% --- Executes during object creation, after setting all properties.
function menuSolver_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menuSolver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in buttonSolve.
function buttonSolve_Callback(hObject, eventdata, handles)
% hObject    handle to buttonSolve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
choice = get(handles.menuSolver,'Value');
global opts;
l_1_str = get(handles.editL1,'String');
l_1 = eval(l_1_str);
l_2_str = get(handles.editL2,'String');
l_2 = eval(l_2_str);
k_0_str = get(handles.editK0,'String');
k_1_str = get(handles.editK1,'String');
k_0 = str2num(k_0_str);
k_1 = str2num(k_1_str);
directions_count_str = get(handles.editDirectionCount,'String');
directions_count_str
directions_count = str2num(directions_count_str);
%opts = struct('ellipses','','approximations','','centers','','directions','');
if (choice == 1) 
    [opts.approximations, opts.centers, opts.ellipses] = ...
        no_good_curves(l_1,l_2,k_0,k_1,directions_count);
elseif  (choice == 2)
    [opts.approximations, opts.centers, opts.ellipses,opts.directions] = ...
        good_curves(l_1,l_2,k_0,k_1,directions_count);
elseif  (choice == 3)
    [opts.approximations, opts.centers, opts.ellipses] = toolbox(l_1,l_2,k_0,k_1,directions_count);
end


% colormap copper

if (get(handles.checkboxHold,'Value') == 1)
    axes(handles.axesApproximation);
    hold(handles.axesApproximation,'on');
    axes(handles.axesEllipses);
    hold(handles.axesEllipses,'on');
else
    axes(handles.axesApproximation);
    hold(handles.axesApproximation,'off');
    axes(handles.axesEllipses);
    hold(handles.axesEllipses,'off'); 
end

if  (choice == 3)
    axes(handles.axesApproximation);
    grid(handles.axesApproximation,'on');
   rotate3d(handles.axesApproximation,'on');
    plot_ia(opts.approximations,'r',struct('fill',0,'shade',0.7));
else
for k = k_0:k_1
   k_shifted = k-k_0+1;
   len = size(opts.approximations{k_shifted}, 2);
   leftk=k-1/2;
   rightk=k+1/2;
   if (k==k_0)
       leftk=k;
   elseif (k==k_1)
       rightk=k;
   end
      [leftk,rightk] 
   axes(handles.axesApproximation);
   colormap copper
   mesh(handles.axesApproximation,...
        linspace(leftk, rightk, len)'*ones(1, len), ...
        ones(1, len)'*opts.approximations{k_shifted}(1, :),...
        ones(1, len)'*opts.approximations{k_shifted}(2, :));
   hold on;
   grid(handles.axesApproximation,'on');
   rotate3d(handles.axesApproximation,'on');
end
 %xlabel('$t$','interpreter', 'latex');
 %ylabel('$l_1$','interpreter', 'latex');
 %zlabel('$l_2$','interpreter', 'latex');
end

axes(handles.axesEllipses);
if  (choice == 3)
    color = 'r'
else 
    color = 'b'
end
for (i = 1:directions_count)
    coords = ellipsoidalProjection(opts.centers{k_1-k_0+1}, ...
        opts.ellipses{k_1-k_0+1,i}, l_1, l_2, 100);
    plot(coords(1,:),coords(2,:),color);
    hold on;
end
%if (handles.menuSolver)

function editK0_Callback(hObject, eventdata, handles)
% hObject    handle to editK0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editK0 as text
%        str2double(get(hObject,'String')) returns contents of editK0 as a double


% --- Executes during object creation, after setting all properties.
function editK0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editK0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editK1_Callback(hObject, eventdata, handles)
% hObject    handle to editK1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editK1 as text
%        str2double(get(hObject,'String')) returns contents of editK1 as a double


% --- Executes during object creation, after setting all properties.
function editK1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editK1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editL1_Callback(hObject, eventdata, handles)
% hObject    handle to editL2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editL2 as text
%        str2double(get(hObject,'String')) returns contents of editL2 as a double


% --- Executes during object creation, after setting all properties.
function editL1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editL2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editL2_Callback(hObject, eventdata, handles)
% hObject    handle to editL2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editL2 as text
%        str2double(get(hObject,'String')) returns contents of editL2 as a double


% --- Executes during object creation, after setting all properties.
function editL2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editL2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editDirectionCount_Callback(hObject, eventdata, handles)
% hObject    handle to editDirectionCount (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDirectionCount as text
%        str2double(get(hObject,'String')) returns contents of editDirectionCount as a double


% --- Executes during object creation, after setting all properties.
function editDirectionCount_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDirectionCount (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxHold.
function checkboxHold_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxHold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxHold
