%> @file DataBaseXls.m
%> @brief 
%> @section matlabComments Details
%> @authors Eoin O'Keeffe (eoin.okeeffe.09@ucl.ac.uk)
%> @date initiated: 13/8/10
%> @date version 1.1: 09/05/2011
%> @version 
%> 1.1
%> @section intro Method
%> Class inherits from DataBaseAbs abstract class.
%> db = DataBaseXls();
%> db.root = '\\VBOXSVR\Dropbox\WP1\WP1 working folder\master model\data';
%> db.db = 'spreadsheetname.xls';
%> db.getAll('worksheetname')
%> @subsection version_history
%> 1.0 
%> <br>1.1 getAll function extended to allow start row number be defined.
%
%> @attention 
%> @todo getAll - add extra parameter (optional) to allow number of columns
%> returned be allowed.
%> @todo getAll - extend so that output can be numeric or just the text
%> data. Not sure how usefule this is. Just an optional output.
classdef DataBaseXls < DataBaseAbs
    %DATABASEXLS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties

    end %properties
    
    methods 
        function open()
            
        end
        function close()
            
        end
        function [num,txt] = rows(table)
            
        end
        function [num,txt] = columns()
        end
        % ======================================================================
    %> @brief getAll read in all data from a worksheet. 
    %>
    %> @param obj instance of the DataBaseXls class.
    %> @param table name of worksheet from which to get the data. Can also
    %> use numbers.
    %> @param startRow Index of row number from which to start (not
    %> required). Useful for files with no column headers
    %> @retval cell_data returned rows from excel/csv file as a cell array
    % ======================================================================
        function [cell_data] = getAll(obj,table,startRow)
            if nargin ==2 
                   startRow = 2; % default to 2 to exclude column headers
            end %if
            %read in rows according to row ids passed
            [numData,txtData,raw] = xlsread(strcat(obj.root,obj.db),table);
            cell_data = raw(startRow:end,:); % Strip off headings
        end %function getAll
        
        function [cell_data] = getColHeaders(obj,table)
            %read in rows according to row ids passed
            [numData,txtData,raw] = xlsread(strcat(obj.root,obj.db),table);
            cell_data = raw(1,:); 
        end %function getColHeaders
        function getCol(table,x)
            
        end %
        function getRow(table,y) 
            
        end %function GetRow
        function [val] = getOne(obj,table,id,col_index)
            [numData,txtData,raw] = xlsread(strcat(obj.root,obj.db),table)
            tmp = cell2mat(raw(2:end,1));
            val = raw(1+(tmp(:)==id),col_index);
            
        end %function 
        function setAll(array)
            
        end %function setAll
        function setCol(x,array)
            
        end %function setCol
        function setRow(y,array)
        end %setROw
        function setOne(x,y, value)
            
        end %setOne
        function saveXls(obj,filename,data,sheet)
            if nargin==3
               sheet = 'sheet2'; 
            end
            xlswrite(strcat(obj.root,filename), data,sheet);
        end
    end %methods
    
end %class DataClass

