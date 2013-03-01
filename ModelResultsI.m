%> @file ModelResultsI.m
%
%> @brief Class inserts results into the ModelResults database
%
%> @section matlabComments Details
%> @authors Eoin O Keeffe (eoin.okeeffe.09@ucl.ac.uk)
%> @date initiated: 21/06/2011
%> <br /> Version 1.1: 22/07/2011
%
%> @version 
%> 1.0:  just added the WriteResults function
%> <br /> Version 1.1: WriteResults title changed to writeResults. Added
%> getResults function that returns results based on an inputted modelrunid.

%
%> @section intro Method
%> 
%
%> @attention 
%> @todo Allow user not to pass ModelRunID, so the funciton just
%> increments the modelrunid from max in the database. ALthough this is
%> probably not useful, as more likely than not, the caller will be
%> inserting multiple rows for a single run.
%> <br />Allow the user create a new model
classdef ModelResultsI < handle
     properties (Constant)
        
     end
      properties 
          end

    methods  (Access=public, Static=true)
        % ======================================================================
    %> @brief createNewModel. Inserts a new model into the database, but
    %also returns the id of an existing model
    %
    %> @param name
    %> @param existing (optional) binary to indicate whether new or
    %> existing model default is new
    %> @param description (optional) for inserting a new model, adding a
    %> description is optional
    %> @param version (optional) for inserting a new model, adding a
    %> version number is optional
    %> @retval modelId
    %> @retval modelRunId
    function createNewModel()
        db = DataBasePG;
        db.db = 'ModelResults';
        
    end %createNewModel
    
% ======================================================================
    %> @brief createNewModel. Inserts a new model into the database, but
    %also returns the id of an existing model
    %
    %> @param name
    %> @param existing (optional) binary to indicate whether new or
    %> existing model default is new
    %> @param description (optional) for inserting a new model, adding a
    %> description is optional
    %> @param version (optional) for inserting a new model, adding a
    %> version number is optional
    %> @retval modelId
    %> @retval modelRunId
    function [modelId,modelRunId] = createNewModelRun(name,existing,description,version)
        if nargin == 1
            existing = 0;
            description = '';
        elseif nargin ==2
            description = '';
        end %if
        db = DataBasePG;
        db.db = 'ModelResults';
        if existing == 0
            % = db.setRow({name 1 date},{'"ID"' '"Name"' '"Version"' '"Description"'},'"Models"');
        else 
            modelId =db.executeQuery(['SELECT "ID" from "Models" where "Name" = ''' name '''']);
        end %if
        % Now create new run based on this
        modelRunId = db.setRow({name modelId date},{'"ID"' '"Description"' '"ModelID"' '"RunDate"'},'"ModelRuns"');
    end %function createNewModel
        % ======================================================================
    %> @brief writeResults - does what it says, writes the results to the
    %> database in the correct format
    %> Function relies on DataBasePG class,Useful class
    %>
    %> @param actualFields - a single row vector cell of the names of the
    %> fields that correspond to the cols in resultsArray. This doesn't have
    %>to be the same number of associated fields in the Models database. In
    %> other words, don't have to have a value for every field.
    %> @param resultsArray cell arry of results to be inputted to database
    %> - can be int, double or string 
    %> @param modelRunID - the ID of the run that results will be
    %> associated with. This should have been generated before the modelling
    %> process began - can be gotten by calling createNewModel above
    %> @retval
    % ======================================================================
    function writeResults(actualFields,resultsArray,modelRunID)
        %First get model details from database
        db = DataBasePG;
        db.db = 'ModelResults';
        %get field names
        fields = db.executeQuery(['Select "Models"."Fields" from "Models"'...
            ' left outer join "ModelRuns" on ("Models"."ID" = '...
            '"ModelRuns"."ModelID") where "ModelRuns"."ID"= '...
            num2str(modelRunID)]);
        %split fields so we can get the insertion column headings
        fields = textscan(char(fields),'%s','delimiter','|')';
        %loop through results and insrt to database
        repFields = cell(size(actualFields,1),1);
        %Build insertion string
        for i=1:size(repFields,1)
            repFields(i,1) = {['"' char(fields{1,1}(Useful.FindCellInCell(fields{1,1},...
                strrep(actualFields(i,1),'"',''))+1,1)) '"']};
        end
        db.setRow(num2cell([repmat(modelRunID,size(resultsArray,1),1) resultsArray]),[{'"ID"' '"ModelRunID"'} repFields'],'"ModelResults"');
        %write to database
    end %function getStorage
           % ======================================================================
    %> @brief GetResults based on a passed modelRunID, get all the results,
    %> resolve them into their correct format and return them to the
    %> caller. A problem with function is that the user needs to know what
    %> order they are gettign stuff in to some extent, although i do pass
    %> back the fields in another cell array
    %>
    %> @param resultsArray cell arry of results to be inputted to database
    %> - can be int, double or string 
    %> @param modelRunID
    %> @param resultName Name of the result series to use, as there can be
    %> multiple result types in fields (optional)
    %> @retval
    % ======================================================================
    function [res actualFields] = getResults(modelRunID,resultName)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Edited EO'K 27/07/2011
        % Added extra input for resultName
        if strcmp(resultName,'')
            resultName = '';
        end %if
        %First get model details from database
        db = DataBasePG;
        db.db = 'ModelResults';
        %get field names
        fields = db.executeQuery(['Select "Models"."Fields" from "Models"'...
            ' left outer join "ModelRuns" on ("Models"."ID" = '...
            '"ModelRuns"."ModelID") where "ModelRuns"."ID"= '...
            num2str(modelRunID)]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Edited EO'K 27/07/2011
        % if resultname not '' then only use fields for this resultname
        if ~strcmp(resultName,'')
            % the results types are split by ||
            fields = textscan(char(fields),'%s','delimiter','||')';
            % Now find the one that corresponds to resultName
            % The one corresponding to resultName is the field following
            % resultName
            for i=1:size(fields{1},1)
               if strcmp(fields{1}(i,1) ,resultName)
                   fields = textscan(char(fields{1}(i+1,1)),'%s','delimiter','|')';
               end %if
            end %for
        else
            %split fields so we can get the insertion column headings
            fields = textscan(char(fields),'%s','delimiter','|')';
        end %if
        %%%%%%%%%%%%%%%%%%%%%
        
        actualFields =  fields{1}(1:2:size(fields{1},1),1);
        repFields =  fields{1}(2:2:size(fields{1},1),1);
        queryStr = 'Select ';
        for i=1:size(repFields,1)
            queryStr = [queryStr '"' char(repFields(i,1)) '",'];
        end %for i
        queryStr =[substr(queryStr,0,size(queryStr,2)-1) ' from "ModelResults" '...
            ' where "ModelRunID" = ' num2str(modelRunID)];
        res = db.executeQuery(queryStr);
        
            
        end%function getResults
    end %methods static
end %classdef

