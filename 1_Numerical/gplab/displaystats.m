function displaystats(message,data,generation)
%DISPLAYSTATS    Display statistical information regarding a GPLAB run.
%   DISPLAYSTATS(MESSAGE,DATA,GENERATION) displays MESSAGE either on the text terminal or,
%   if open, on the Graphical User Interface of GPLAB. DATA and GENERATION are variables
%   used in MESSAGE.
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
       message2gui = sprintf(message,data);
       message2gui = strtrim(message2gui);
                    
       if gui_handle.Generation ~= generation
           infomessage = gui_handle.tempMessages{1};
           for i = 2:length(gui_handle.tempMessages);
               infomessage = [infomessage sprintf(['\n' gui_handle.tempMessages{i}])];
           end
           set(output_handle.statisticsInfo,'visible','on','String', infomessage );
           
           gui_handle.messageCounter = 1;
           gui_handle.tempMessages = {}; 
       end
  
       gui_handle.tempMessages{gui_handle.messageCounter} = message2gui;
       gui_handle.messageCounter = gui_handle.messageCounter + 1;
       
       drawnow
       gui_handle.Generation = generation;
       guidata(gui_figure,gui_handle);  
   else
       fprintf(message,data);
   end
   
end

