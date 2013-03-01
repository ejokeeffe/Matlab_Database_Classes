%> @file DataBaseAccess.m
%> @brief Abstract class for inheritance for DataBase subclasses (eg for
%xls, access, mysql etc)

%> @mainpage Overview

%> @section intro Introduction

%> Authors: Eoin O'Keeffe,

%> Date initiated: 13/8/10

%> Version: 1.0

%> This class provides access to database files of type xls/xlsx
classdef DataBaseAccess < DataBaseAbs
    %DATABASEACCESS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        con
    end %properties
    
    methods 
        function [success] = open(obj)
            % Specify the location of the database on disk.
                dbpath = strcat(obj.root,obj.db);

            % Create the connection URL.
            conurl = ['jdbc:odbc:Driver={Microsoft Access Driver (*.mdb)};' ...
            'DBQ=' dbpath];

            % Connect to the database.
            obj.con = database('','','','sun.jdbc.odbc.JdbcOdbcDriver', conurl);
    
        end
        function close()
            
        end
        function [num,txt] = rows(table)
            
        end
        function [num,txt] = columns()
        end
        
        function [cell_data] = getAll(obj,table) 
            obj.open;
            e = exec(obj.con, ['select * from ' table]);
            %this should return as cell array, if not then change in 
            %setdbprefs to DataReturnFOrmat : 'cellarray'
            cell_data = fetch(e);
            cell_data = cell_data.Data;
            close(obj.con);
        end %function getAll
        function [cell_data] = executeQuery(obj,sqlstr);
            obj.open;
            e = exec(obj.con,sqlstr);
            cell_data=fetch(e);
            cell_data = cell_data.Data;
            close(obj.con);
        end %function
        %where a table links to another table 
        function [cell_data] = getLinkedData(obj,parenttable,childtable,parentlink,childlink)
            
        end %function getLinkedData
        function getCol(table,x)
            
        end %
        function [ID] = getID(obj,sqlstr,table)
            obj.open;
            e = exec(obj.con, ['select ID from ' table ' ' sqlstr]);
            %this should return as cell array, if not then change in 
            %setdbprefs to DataReturnFOrmat : 'cellarray'
            cell_data = fetch(e); 
            close(obj.con);
            ID = cell_data.data(1,1);
            
                if strcmp(ID,'No Data')
                    ID = -1;
                else
                    ID = cell2mat(ID);
                
                end %If
            
            
        end %function getID
        function getRow(table,y) 
            
        end %function GetRow
        function [val] = getOne(obj,table,id,col_index)
            
            
        end %function 
        function setAll(array)
            
        end %function setAll
        function setCol(obj,array,colnames,table,ID)
            obj.open();
            update(obj.com,table,colnames,array,['WHERE ID =' ID]);
            close(obj.con);
        end %function setCol
        function update(obj,array,colnames,table,ID)
            obj.open();
            update(obj.con,table,colnames,array,['WHERE ID =' Useful.Val2Str(ID)]);
            close(obj.con);
        end %function Update
        function [ID] = setRow(obj,array,colnames,table)
            %con = database(obj.db, '', '');
            obj.open();
            %TODO This should be improved but get the max ID and then increment
            %by 1 to insert
            e = exec(obj.con,['SELECT MAX(ID) FROM ' table]);
            ID = fetch(e); 
            ID=cell2mat(ID.data) + 1; %Returns a structure, the data is contained within the .data variable
            if isnan(ID)
                ID = 1;
            end %if
            %Not build sql string
            %tmp = '';
            %for i = 1:size(array,2)
            %    tmp = strcat(tmp,',''',array(i),'''');
            %end %for
            %tmp = strcat(ID(1,1),tmp);
            %table = [' ' table ' ']; %have to add spaces, but strcat trims it
            %TODO tidy this up, matlab is really crap at handling strings
            %sql_str = strcat('insert into ',table,' VALUES (',tmp,')');
            %run sql statement
            %res = exec(obj.con, sql_str);
            
            insert(obj.con,table,colnames,[ID array]);
            close(obj.con);
        end %setROw
        function setOne(x,y, value)
            
        end %setOne
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> @brief get the column names from a atable
        %> @param obj
        %> @param tbleName
        %> @param addQuotes
        %> @retval 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [cols] = getColumnNames(obj,tbleName,addQuotes)
            sql = sprintf(['select column_name '...
                'from information_schema.columns '...
                'where table_name = ''%s'''],tbleName);
            cols = flipud(obj.executeQuery(sql))';
            if nargin>2
                if addQuotes==1
                    for i=1:numel(cols)
                        cols(i) = {['"' char(cols(i)) '"']};
                    end %for
                end %if
            end %if
        end %function getCOlumnNames
    end %methods
    
end %class DataClass

