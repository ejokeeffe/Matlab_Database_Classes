%> @file DataBasePost.m
%> @brief Class allows access to a PostGres server
%> @section matlabComments Details
%> @authors Eoin O'Keeffe (eoin.okeeffe.09@ucl.ac.uk)
%> @date initiated: 08/03/2011
%> @date 1.1: 23/05/2011
%> @date 1.2: 17/06/2011
%> @data 1.3: 22/07/2011
%> <br /> Version 1.4: 16/08/2011
%> <br /> <i>Version 1.5</i>: 23/04/2012
%> <br /> <i>Version 1.6</i>: 6/11/2012
%> <br /> <i>Version 1.7</i>: 08/01/2012
%> <br /> <i>Version 1.8</i>: 05/03/2013
%>
%> @version 
%> 1.1 New class created as matlab doesn't connect with Access in the
%> 64-bit version. Handy that isn't it.
%><br /> Version 1.1: AddRow edited so an array of rows can be passed and
%> each is inserted in the database
%><br /> Version 1.2: Edited setRow so there is no id field and a straigh
%>insert is done
%><br /> Version 1.3: extended getAll function so it has an extra optional
%> variable, so the user can append a querysting on the end of the get all
%> call
%><br /> Version 1.4: executeQuery adjusted so it bruns multiple queries if
%> the number of rows returned is very large
%> <br /> <i>Version 1.5</i>: in getAll function, all users specify which
%> cols to return
%> <br /> <i>Version 1.6</i>: Added runBlockInserts to batch insert a group
%> of inserts at one time rather than individually. In time, setRow should
%> be amended to incorporate this function and then make this function.
%> Orderby field also added to executeQuery, although this was done some time
%> ago
%> <br /><i>Version 1.7</i>: added function to escape a string value for
%> insert into the database
%> <br /><i>Version 1.8</i>: Add getting useranme and password from config
%> file

%> private
%> @section intro Method
%> Usual database methods
%
%> @attention For the insert statement i have to loop through the datarows,
%> insert won't allow be insert multiple rows at once
%> @todo setRow - Add code so not to increment the ID
classdef DataBasePG < DataBaseAbs
    %DATABASEACCESS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        con
        pass
        user
    end %properties
    
    methods 
        function obj = DataBasePG()
            %check if java on classpath
            staticPath = javaclasspath('-static');
            dynamicPath = javaclasspath('-dynamic');
            if 0 == (sum(cellfun(@(x) ~isempty(regexp(x,'postgres')),staticPath)) + ...
                    sum(cellfun(@(x) ~isempty(regexp(x,'postgres')),dynamicPath)))
                    fname = mfilename('fullpath');
                    fname = fname(1:regexp(fname,mfilename)-1);
                    javaaddpath(fullfile(fname,'postgresql-9.1-901.jdbc4.jar'));
            end %if
        end %constructor
        function [success] = open(obj)
            % Create the connection URL.
            conurl = ['jdbc:postgresql://localhost:5432/' obj.db];

            % before connecting to database, make sure we have a password
            % and username
            %retrieve from config file
            obj.user = Useful.getConfigProperty('dbuser');
            obj.pass = Useful.getConfigProperty('dbpassword');
            
            % Connect to the database.
            obj.con = database(obj.db,obj.user,obj.pass,'org.postgresql.Driver', conurl);
    
        end
        function close()
            
        end
        function [num,txt] = rows(table)
            
        end
        function [num,txt] = columns()
        end
        % ======================================================================
    %> @brief gets all rows from a table, can also add a querystring on the
    %> end as an option. Edited EOK 23rd April 2012 to allow user pass the
    %> cols to get from the database
    %>
    %> @param obj
    %> @param table table name not including "" if it has capitals
    %> @param queryString (optional) in the form 'where "X" = 10' - char
    %> format obviously
    %> @param cols - string array of cols for the get all statement
    %> @retval cell_data all returned data from the database query as a
    %> cell array
    % =====================================================================
        function [cell_data] = getAll(obj,table,queryString,cols)
            if nargin == 2
                queryString = '';
                cols = '*';
            elseif nargin == 3
                cols = '*';
            else
                tmp = '';
                for i=1:numel(cols)
                    tmp =[tmp ',' cols{i}];
                end %for i
                cols = tmp(2:end);
            end %if
            obj.open;
            e = exec(obj.con, ['select ' cols ' from "' table '" ' queryString]);
            %this should return as cell array, if not then change in 
            %setdbprefs to DataReturnFOrmat : 'cellarray'
            cell_data = fetch(e);
            cell_data = cell_data.Data;
            close(obj.con);
        end %function getAll
        % ======================================================================
    %> @brief executeQuery - execute a sqlstr. Also allows the query to be
    %> split
    %>
    %> @param obj instance of the DataBasePG class.
    %> @param sqlstr
    %> @param countQuery the sqlstr for count the rows returned (optional)
    %> @param orderByField required so the offsetting and limit sqlquery
    %> works properly
    %> @retval cell_Data 
    % ======================================================================
        function [cell_data] = executeQuery(obj,sqlstr,countQuery,orderByField)
            
            if nargin==2
                %Just a normal fetch
               obj.open;
                e = exec(obj.con,sqlstr);
               cell_data=fetch(e);   
                cell_data = cell_data.Data;
                close(obj.con);
            else
                %we may have a lot of rows to return, set the limit to
                %30,000
                obj.open;
                e = exec(obj.con,countQuery);
                countRows=fetch(e);
                countRows = cell2mat(countRows.Data);
                %now loop through and add in results until we hit the limit
                exitWhile = 0;
                limit = 30000;
                offset = 0;
                cell_data = [];
                while exitWhile == 0
                    e = exec(obj.con,[sqlstr '  ORDER by ' orderByField ' LIMIT ' num2str(limit) ...
                        ' OFFSET ' num2str(offset)]);
                    tmpCell = fetch(e);
                    tmpCell = tmpCell.Data;
                    cell_data = [cell_data;tmpCell];
                    %increment offset
                    offset = limit +offset;
                    if countRows< limit+offset
                        if countRows<=offset
                            exitWhile=1;
                        else
                            limit = countRows - offset;
                        end %if
                    end %if     
                end %while
                close(obj.con);
            end %if
            
        end %function
        %where a table links to another table 
        function [cell_data] = getLinkedData(obj,parenttable,childtable,parentlink,childlink)
            
        end %function getLinkedData
        function getCol(table,x)
            
        end %
        function [ID] = getID(obj,sqlstr,table)
            obj.open;
            e = exec(obj.con, ['select ID from "' table '" ' sqlstr]);
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
        % ======================================================================
    %> @brief inserts a new row into a database table
    %>
    %> @param obj
    %> @param array array of values to insert that correspond to the
    %> colnames
    %> @param colnames vector of colnames (cell array obviously)
    %> @param table name of table to insert into
    %> @param noID optional parameter to be set to 0 if there is no field
    %> in the table that is called ID, so the row is inserted without
    %> incrementing an ID field
    %> @param fastInsert - uses fastinsert instead of insert
    %> @retval ID of row just inserted
    % =====================================================================
        function [ID] = setRow(obj,array,colnames,table,noID,fastInsert)
            if nargin == 5
               %Add code to not increment the ID as there is no ID field there 
               if noID == 1
                   %Set ID to zero, wo when it is decremented at the end it
                   %returns -1.
                   ID=0;
               else 
                   noID =0;
               end %if
               fastInsert = 0;
            elseif nargin ==6
                %Add code to not increment the ID as there is no ID field there 
               if noID == 1
                   %Set ID to zero, wo when it is decremented at the end it
                   %returns -1.
                   ID=0;
               else 
                   noID =0;
               end %if
            else
                fastInsert = 0;
                noID = 0;
            end
            %increment the ID
            %con = database(obj.db, '', '');
            obj.open();
            %TODO This should be improved but get the max ID and then increment
            %by 1 to insert
            %Edited EO'k 17/06/2011 so only do ID thing if there is an ID
            %field
            if noID==0
                e = exec(obj.con,['SELECT MAX("ID") FROM ' table]);
                ID = fetch(e); 
                ID=str2num(Useful.Val2Str(ID.data)) + 1; %Returns a structure, the data is contained within the .data variable
                if isnan(ID)
                    ID = 1;
                end %if
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
            %Added EOK 23/05/2011, add row for each row in array
            %Edited EOK 15/04/2011, to allow ability to add all entries at
            %once
            if fastInsert == 0
                for i=1:size(array,1)
                    i;
                    if noID==0
                        insert(obj.con,table,colnames,[{ID} array(i,:)]);
                        ID= ID+1;
                    else 
                        insert(obj.con,table,colnames,[array(i,:)]);
                    end
                end
            else
                %use fastinsert instead
                for i=1:size(array,1)
                    i
                    if noID==0
                        fastinsert(obj.con,table,colnames,[{ID} array(i,:)]);
                        ID= ID+1;
                    else 
                        fastinsert(obj.con,table,colnames,[array(i,:)]);
                    end
                end
            end
            %Step back to last id set
            ID=ID-1;
            close(obj.con);
        end %setROw
        function setOne(x,y, value)
            
        end %setOne
        % ======================================================================
    %> @brief gets the column titles from a given datatable
    %>
    %> @param obj
    %> @param tbleName tablename
    %> @param addQuotes 1 or 0 or empty, if 1 then double quotes are added
    %> around col titles
    %>
    %> @retval returns a cell array of the column titles
    % =====================================================================
        function [cols] = getColumnNames(obj,tbleName,addQuotes)
            sql = sprintf(['select column_name '...
                'from information_schema.columns '...
                'where table_name = ''%s'' order by ordinal_position'],tbleName);
            cols = obj.executeQuery(sql);
            if nargin>2
                if addQuotes==1
                    for i=1:numel(cols)
                        cols(i) = {['"' char(cols(i)) '"']};
                    end %for
                end %if
            end %if
        end %function getCOlumnNames
        
        % ======================================================================
    %> @brief gets the column titles from a given datatable
    %>
    %> @param obj
    %> @param tbleName tablename
    %> @param addQuotes 1 or 0 or empty, if 1 then double quotes are added
    %> around col titles
    %>
    %> @retval returns a cell array of the column titles
    % =====================================================================
        function [tables] = getTableNames(obj)
            sql = ['select distinct(table_name) from information_schema.columns'...
                ' where table_schema = ''public'' order by table_name'];
            tables = obj.executeQuery(sql);
            
        end %function getTableNames
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> @brief this is additional to setRow, it used to run a batch of
        %> insert statements as one. We don't need to run them individually
        %>
        %> @param obj
        %> @param fields a cell array of characters with the fieldnames 
        %> @param table name of table exlcluding double quotes
        %> @param vals matrix of values to insert
        %> @param insertCount [optional] can set the number of rows to
        %> insert each execute. Defaults to 200
        function [success] = runBlockInserts(obj,fields,table,vals,insertCount)
            %set the inserCount
            if nargin==4
                insertCount=200;
            end %if
            %build the fields string
            sqlFields = '';
            for i=1:length(fields)
                sqlFields = sprintf('%s,%s',sqlFields,fields{i});
            end %for i
            sqlFields = sqlFields(2:end);
            %Now build the values
            %and run insert
            insertStatement ='';
            for i=1:size(vals,1)
                %build the values string
                sqlValues = '';
                for j=1:size(vals,2)
                    if isa(vals{i,j},'char')
                        sqlValues = sprintf('%s,''%s''',sqlValues,vals{i,j});
                    else
                        sqlValues = sprintf('%s,%d',sqlValues,vals{i,j});
                    end %if
                    
                end %for j
                sqlValues = sqlValues(2:end);
                %build the sql statement    
                insertStatement = sprintf(['%s Insert into "%s" (%s) values'...
                    ' (%s);'],insertStatement,table,sqlFields,sqlValues);
                %do we run it now or what?
                if (rem(i,insertCount)==0 || i==size(vals,1))
                    %run the statement
                    obj.executeQuery(insertStatement);
                    insertStatement='';
                end %if
            end %for i
        end %runBlockInserts
        
    end %methods
    
    methods(Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> @brief escapes restricted operators for a postgres sql statement
        %>
        %> @param value the value to escape. It's check if it's a string
        %>
        %> @retval the escape value returned
        function [escaped] = escapeValue(value)
            escaped =value;
            if isa(value,'char')
                escaped = strrep(value,'''','''''');
            end %char
        end %function escapeValue
    end %static methods
end %class DataClass

