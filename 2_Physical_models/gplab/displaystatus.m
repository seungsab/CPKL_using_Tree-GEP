function displaystatus(message,varargin)
%DISPLAYSTATUS    Display informations and warnings regarding a GPLAB run.
%   DISPLAYSTATUS(MESSAGE) displays MESSAGE either on the text terminal or,
%   if open, on the Graphical User Interface of GPLAB.
%
%   Created (2014) by Nuno Tenazinha
%   This file is part of the GPLAB Toolbox

   %gui_figure = findobj('Name','GPLAB Graphical User Interface');
   if exist('gui_figure','var') && ~isempty(gui_figure)
       gui_handle = guidata(gui_figure);
       output2find = ['Output Messages ' '(' num2str(gui_handle.runnumber(end)) ')'];
       output_figure = findobj('Name', output2find);
   end
   
   if exist('output_figure','var') && ~isempty(output_figure)
       output_handle = guidata(output_figure); 

       message = strrep(message,'\n','');
       message = strrep(message,'- ','');
       message = strrep(message,' -','');
       message2gui = sprintf(message,varargin{:});
       
       if(nargin > 1 && ~isempty(strfind(message, 'generation')))
           if isempty(gui_handle.loadedRun)
               progress = varargin{1}/gui_handle.gnumber;
           else
               progress = (varargin{1}-gui_handle.loadedRun.vars.state.generation)/gui_handle.gnumber;
           end
           set(output_handle.progressBar,'Position',[0 0.1 progress 0.80],'String',[num2str(round(progress*100)) '%']);
       end
       output_handle.outputMessagesList{end+1} = message2gui;    
       messagenumber = length(output_handle.outputMessagesList);
       set(output_handle.outputMessages, 'String', output_handle.outputMessagesList)% flipud(output_handle.outputMessagesList')')
       set(output_handle.outputMessages,'ListboxTop',messagenumber);
       guidata(output_figure,output_handle) 
       drawnow
   else
       fprintf(message,varargin{:});
   end

   clear output_figure gui_figure output_handle gui_handle scale
end