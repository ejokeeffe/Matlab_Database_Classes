%> @file Useful.m
%> @brief Class built to do parsing etc containing handy functions
%> @section matlabComments Details
%> @authors Eoin O'Keeffe (eoin.okeeffe.09@ucl.ac.uk)
%> @date initiated: 10/10/2010
%> <br />version 1.1: 01/05/2011
%> <br />version 1.2: 13/05/2011
%> <br />version 1.3: 20/07/2011
%> <br /><i>version 1.4:</i> 10/11/2011
%> <br /><i>version 1.5:</i> 12/03/2011
%> <br /><i>version 1.6:</i> 20/01/2013
%> <br /><i>version 1.7</i>: 06/02/2013
%>
%> @version 
%> 1.0 
%> <br>1.1 Added conversion to upper case for FindCellInCellSimple
%> <br />Version 1.2 Added function to build matrix for viewing in excel
%> from a database input. It's a reverse versions of convertsquareto3col
%> <br /> Version 1.3: Added matchUniqueArrays function
%> <br /><i> Version 1.4:</i> Added convertNColToMatrix function that
%> converts a n column 2 d array into a n-1 dimensional array. This should
%> supercede convert3colToMatrix, so that function is now popping out a
%> warning
%> <br /><i> Version 1.5:</i> Added function to parse lat and lon string
%> variables to return two double with minus for North or west
%> <br /><i> Version 1.7</i>: Addition of direction param in
%> GreatCircleDistance to get the distance in a particular direction
%>
%> @section intro Method
%> All methods are static as they're effectively general parsing function.
%> It's effectively a repositiory for functions. It's not a real class.
%> <br /><i> Version 1.6:</i> Addition of a function that converts a value
%> to a char for writing to an excel file (ie getting the column field)
%> @subsection version_history
%
%> @attention 
classdef Useful
    %USEFUL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
    end
    methods (Static = true) % Define static method 
     function [x] = Val2Str(val)
         if iscell(val)
            if cellfun('isclass',val,'double') | cellfun('isclass',val,'single')
                x = num2str(cell2mat(val));
            else
                x=char(val);
            end %if
         elseif isnumeric(val)
             x=num2str(val);
         else
             %just set as input value
             x=char(val);
         end %if
     end % grid
     function [x] = Cell2Str(cell)
         x = mat2str(cell2mat(cell));
     end %function Cell2Str
     function [x] = FindCellInCell(searchCell,elementCell,indexes,fullmatch)
         %The order of indexes should be the same order as searchCell, but
         %its cell values should correspond to the cell in elementCell
         if nargin >= 3
            %the indexes have been passed, these should be in a cell array
            %but they may not all be values
            tmp_elementcell = 0;
            tmp_searchcell = 0;
            for i=1:size(indexes,2)
                if ~strcmp(Useful.Val2Str(indexes(1,i)),'[]') &...
                        ~strcmp(Useful.Val2Str(indexes(1,i)),'-1')
                    %add index to tmp_cell
                    tmp_elementcell = [tmp_elementcell str2double(cell2mat(indexes(1,i)))];
                    tmp_searchcell = [tmp_searchcell i];
                end %if
            end %for
            %now strip off the unwanted indexes in element cell
            elementCell = elementCell(1,tmp_elementcell(1,2:end));
            searchCell = searchCell(:,tmp_searchcell(1,2:end));
         end %if
         fuzzymatch = 0;
         if nargin == 4
             if fullmatch > 0 
                 fuzzymatch = 1;
             end %if
         end %if
         %Now lets make sure we're comparing cell array of the same cell
         %classes - ie strings
         searchCell = cellfun(@Useful.Val2Str,searchCell,'Un',0);
         elementCell = cellfun(@Useful.Val2Str,elementCell,'Un',0);
         x = -1;
         if fuzzymatch ==0
             for i=1:size(searchCell,1)
                 if isequal(searchCell(i,:),elementCell)
                    x = i; 
                    break;
                 end
             end %for
         else
             for i=1:size(searchCell,1)
                 %loop through each item in the column fields and check if
                 %there are any matches
                 for j=1:size(searchCell,2)
                     if isequal(searchCell(i,j),elementCell(j))
                        x = i; 
                        break;
                     end %if
                 end %for
             end %for
         end %if
     end %FindCellInCell
     
     function [x] = FindCellInCellSimple(searchCell,elementCell)
         x=0;
                  %Now lets make sure we're comparing cell array of the same cell
         %classes - ie strings
         %EOK amended 5th may 2011 to convert both to upper case to allow
         %comparison
         searchCell = cellfun(@Useful.Val2Str,searchCell,'Un',0);
         searchCell = cellfun(@upper,searchCell,'Un',0);
         elementCell = cellfun(@Useful.Val2Str,elementCell,'Un',0);
         elementCell = cellfun(@upper,elementCell,'Un',0);
         for i=1:size(searchCell,1)
                 if isequal(searchCell(i,:),elementCell)
                    x = i; 
                    break;
                 end
             end %for
     end %FindCellInCellSimple
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %> @brief This not use this function anymore - use the inbuilt
     %> distance function in matlab
     %> @param sourceLon
     %> @param sourceLat
     %> @param destLon
     %> @param destLat
     %> @retval dist_nm 
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     function [dist_km] = GreatCircleDistance(sourceLon, sourceLat, ...
             destLon, destLat,direction)
         if nargin == 4
             direction = [];
         end %if
         %disp(size(phi_f));
         %disp(['x:' num2str(phi_s) ' ' num2str(lambda_s) ' y:' num2str(phi_f) ' ' num2str(lambda_f)]);
         dist_km= zeros(size(destLon,1),1);
         for i=1:size(destLon,1)
%         cos_delta = sin(phi_s*pi/180)*sin(phi_f(i)*pi/180)+...
%             cos(phi_s*pi/180)*cos(phi_f(i)*pi/180)*cos((lambda_s-lambda_f(i))*pi/180);
%         delta_rads = acos(cos_delta);
%changed by eok to formula on movable-type.co.uk/scripts/latlong.html
        r=6371; %kms
            dlat = degtorad((destLat(i)-sourceLat));
            dlon = degtorad((destLon(i)-sourceLon));
            lat1 = degtorad(sourceLat);
            lat2 = degtorad(destLat(i));

            a = sin(dlat/2)^2 + sin(dlon/2)^2*cos(lat1)*cos(lat2);
            c = 2*atan2(sqrt(a),sqrt(1-a));
            if isempty(direction)
                dist=r*c;
            elseif direction == 'E'
                %check if we're going in the wrong direction
                if sourceLon>destLon
                    dist=r*(2*pi-c);
                else
                    dist=r*c;
                end %if
            elseif direction == 'S'
                %check if we're going in the wrong direction
                if sourceLat<destLat
                    dist=r*(2*pi-c);
                else
                    dist=r*c;
                end %if
            elseif direction =='W'
                if sourceLon>destLon
                    dist=r*(2*pi-c);
                else
                    dist=r*c;
                end %if
             elseif direction =='N'
                if sourceLat<destLat
                    dist=r*(2*pi-c);
                else
                    dist=r*c;
                end %if   
            end   
                
            
            %dist_km(i) = dist_km*.539957;
            dist_km(i) = dist;
         
        %disp(['Returning: ' num2str(dist_nm)]);
         end 
         %disp(sprintf('about to exit and inputs are %d and %d. Result is %f',numel(sourceLon),numel(destLon),dist_km(1)));
     end
     %>@brief loops through a lit of lon/lat and gets the gc between them
     function dist = getSequentialGC(lon,lat)
         dist = 0;
         for i=1:numel(lon)-1
             dist = dist + Useful.GreatCircleDistance(lon(i),lat(i),lon(i+1),lat(i+1));
         end %for i
     end %getSequentialGC
     %@brief wrapper function for distance calc
    function [dists] = getMultipleDistances(sourcePoints,destPoints)
        dists = zeros(size(sourcePoints,1)*size(destPoints,1),1);
        dists = Useful.GreatCircleDistance(sourcePoints(:,1), sourcePoints(:,2),...
                         destPoints(:,1), destPoints(:,2));
%         for i=1:size(sourcePoints,1)
%             for j=1:size(destPoints,1)
%                 dists((i-1)*size(destPoints,1) + j,1)=...
%                     distance(sourcePoints(i,2), sourcePoints(i,1),...
%                         destPoints(i,2), destPoints(i,1),...
%                         almanac('earth', 'wgs84'))
%                     
%             end %for j
%         end %for i
    end %getMultipleDistances

%==========================================================================
%> @brief gets the index of the min dist from each point in array 1 to each
%> set of points in array 2
%>
%> @param sourcePoints - dataset with .lon and .lat fields
%> @param destPoints - datset with .lon and .lat fields
%> @param useGreatCircle - indicates whether to use great circl
%> defualt is yes
%> @retval indxsOfMinDist vector of same length as sourcePoints with
%> indexes amthcing destPoints
function [minDist,indxsOfMinDist] = getMinDistBetweenPoints(sourcePoints,destPoints,...
        useGreatCircle)
    if nargin ==2
        useGreatCircle = 1;
    end %if
    %may need to split up the source points - do it in groups of 1000
    if size(sourcePoints,1)>1000
        splitSize = 5000;
        totals = [1:splitSize:size(sourcePoints,1) size(sourcePoints,1)];
        minDist =[];
        indxsOfMinDist =[];
        for i=1:length(totals)-1
            endVal=totals(i+1)-1;
            if i==length(totals)-1
                endVal=totals(i+1);
            end %if
             [dists] = pdist2([destPoints.lon destPoints.lat],...
                 [sourcePoints.lon(totals(i):endVal,:) ...
                 sourcePoints.lat(totals(i):endVal,:)],'euclidean'); 
             %get indx of min
             [tmpMinDist,tmpMinDistIndxs] = min(dists);
             minDist = [minDist;tmpMinDist(:)];
             indxsOfMinDist = [indxsOfMinDist;tmpMinDistIndxs(:)];
        end %for i
      else
              [dists] = pdist2([destPoints.lon destPoints.lat],...
             [sourcePoints.lon sourcePoints.lat],'euclidean'); 
         [minDist,indxsOfMinDist] = min(dists);
    end %if
     indxsOfMinDist = indxsOfMinDist(:);
     minDist=minDist(:);
end %function getMinDistBetweenPoints
     % ======================================================================
    %> @brief COnvert square matrix to three col matrix
    %>
    %> @param sqMat instance of the Commodities class.
    %> @retval x three col matrix in form [origin(1) dest(2) val(3)]
    % ======================================================================
    function [x] = convertSquareTo3Col(sqMat)
        x = cell(size(sqMat,1),3);
        for i=1:size(sqMat,1)
            tmp_origin = repmat(sqMat(i,1),size(sqMat,1)-1,1);
            tmp = sqMat(i,2:end)';
            tmp_dest = sqMat(:,1);
            if i==1 
                tmp = tmp(2:end,:);
                tmp_dest = tmp_dest(2:end,:);
            elseif i==size(sqMat,1)
                tmp = [tmp(1:i,:);tmp(i:end,:)];
                tmp_dest = [tmp_dest(1:i,:);tmp_dest(i:end,:)];
            else
                tmp = tmp(1:end-1,:);
                tmp_dest = tmp_dest(1:end-1,:);
            end
            x(((i-1)*size(sqMat,1))+1:((i-1)*size(sqMat,1))+size(sqMat,1),:)=...
                [tmp_origin tmp tmp_dest];
            
        end
    end
    % ======================================================================
    %> @brief Convert 3 col matrix to square matrix
    %>
    %> @param dta 3 col matrix
    %> @param units value to put in at top left cell
    %> @retval x square matrix
    % ======================================================================
    function [x] = convert3ColToSquare(dta,units)
        if nargin==1
           units=''; 
        end
        %Get unique origin and dests
        %A value in the dest row may not appear in the origin row so append
        %it and get unique
        uniqueVals = unique([dta(:,1);dta(:,2)]);
        x = cell(size(uniqueVals,1),size(uniqueVals,1));
        
        for i=1:size(dta,1)
            disp(['Inserting row ' num2str(i) ' into square matrix']);
            disp(['From ' Useful.Val2Str(dta(i,1)) ' to ' ...
                Useful.Val2Str(dta(i,2))]);
             tmpRow = cellfun(@(x)strcmp(x,dta(i,1)),uniqueVals(:,1));
             tmpCol = cellfun(@(x)strcmp(x,dta(i,2)),uniqueVals(:,1));   
            x(tmpRow,tmpCol)= dta(i,3);
            end
        %Now add column names and row names
        x = [uniqueVals x];
        x = [{units} uniqueVals';x];
    end
    % ======================================================================
    %> @brief Takes in two vectors and rearranges the first so that it's in
    %> line with the second. The first can have more than the second
    %> @param teRearrange vector to rearrange
    %> @param matchArray - vector to rearrange to
    %> @retval indxs - vector that can be used to sort toRearrange
    % =====================================================================
    function [indxs] = matchUniqueArrays(toRearrange,matchArray)
        indxs = zeros(size(matchArray,1),1);
        for i=1:size(matchArray,1)
            matchArray(i,1)
            if isnan(matchArray(i,1)) || matchArray(i,1)==-1 ...
                    || matchArray(i,1)==0
                indxs(i,1) = -1;
            else
                indxs(i,1) = find(toRearrange==matchArray(i,1),1);
            end %if
        end %for
        
    end %function matchUniqueArrays
    % ======================================================================
    %> @brief Takes in two matrices with mathcing no of columns, and
    %> selects the indxs from the first that all match with the second. If
    %> there are multiple matches then the match value is copied for each
    %> item it matches
    %> @param teRearrange vector to rearrange
    %> @param matchMatrix - vector to rearrange to
    %> @retval indxsRearrange - vector that can be used to sort toRearrange
    %> @retval indxsMatch - indxs of match items as there may be duplicates
    % =====================================================================
    function [indxsRearrange indxsMatch] = matchMatrices(toRearrange,matchMatrix)
        indxsMatch = [];
        indxsRearrange = [];
        indxsTemplate = [1:size(toRearrange,1)]';
        for i=1:size(matchMatrix,1)
            tmpRearrange = [];
            for j=1:size(toRearrange,2)
                tmpRearrange = [tmpRearrange toRearrange(:,j)==matchMatrix(i,j)]; 
            end %for j
            % we should have tmpRearrage in the form [1 0 0 1;0 0 0 1] for example,
            % there each value correspondes to whether there is a match for
            % that corresponding col and the total number of rows is the
            % total number of rows in toRearrange
            %Now add our results to indxRearrange
            tmpMatches = ismember(tmpRearrange,ones(size(tmpRearrange,1),size(tmpRearrange,2)),...
                'rows');
            tmpMatches = indxsTemplate(tmpMatches,1);
            indxsRearrange = [indxsRearrange;tmpMatches];
            %How many matchMatrix values to add?
            indxsMatch = [indxsMatch;repmat(i,size(tmpMatches,1),1)];
        end %for i
    end %function matchMatrices
     % ======================================================================
    %> @brief Takes in two variables. The first is the matrix to be
    %aggregated, with multiple rows. 
    %> For each variable in each row in bounds, we sum all the
    %corresponding values in disagMatrix. For example, if we want all the
    %years imports to USA, UK in all the years, bounds will contain a
    %structure with two fields (one of years and one of countries)
    %> @param disagMatrix Matrix containing all the disagregated data
    %> @param bounds unique variables overwhich to sum (structure)
    %> @retval res 
    % =====================================================================
    function [res] = sumOverVariables(disagMatrix,bounds)
        fnamesBounds = fieldnames(bounds);
        %Get total number of fields in res
        rowsRes = [];
        for i=1:numel(fnamesBounds)
            uniqBounds = unique(bounds.(fnamesBounds{i}));
            if isempty(rowsRes)
                rowsRes = uniqBounds;
            else 
                tmpSize = size(rowsRes,1); 
                %Need to create a copy of rowsRes for each of the current 
                %current bounds and reshape 
                tmpRowsRes = [];
                for j=1:size(rowsRes,2)
                    tmpRowsRes = [tmpRowsRes reshape(repmat...
                        (reshape(rowsRes(:,j),1,size(rowsRes,1))...
                        ,size(uniqBounds,1),1),size(rowsRes,1)*size(uniqBounds,1),1)]; 
                end
                %Now append this reshaped rowsRes (ie. tmpRowsRes) to a
                %stacked copy of current uniqBounds
                rowsRes = [tmpRowsRes,repmat(uniqBounds,tmpSize,1)];
            end %if
            
        end %for i
        %Set res in size
        res = zeros(size(rowsRes,1),numel(fnamesBounds)+1);
        %Now drop in summed values
        %Increment last value first
        for i=1:size(res,1)
            % FInd matching rows first
            indxsToSum = Useful.matchMatrices(disagMatrix(:,1:end-1),rowsRes(i,:)); 
            %Now sume over these rows
            res(i,:) = [rowsRes(i,:) sum(disagMatrix(indxsToSum,end))];
        end
    end %function sumOverVariables
    % ======================================================================
    %> @brief converts 3 col to matrix. Eg [category year value] converts
    %> to category value with a column for each of the years. It is assumed
    %> the category is the first col, the second col is the one that'll move
    %> to colums and the third col is the value field
    %>
    %> @param dta 3 col matrix
    %> @retval x matrix
    % ======================================================================
    function [x,categories,fields] = convert3colToMatrix(dta)
       %get unique column fields
       fields = sort(unique(dta(:,2)),'ascend');
       categories = sort(unique(dta(:,1)),'ascend');
       x = zeros(size(categories,1),1+size(fields,1));
       x(:,1) = categories;
       for i=1:size(dta,1)
           indxCategory = find(categories(:,1)==dta(i,1));
           indxField = find(fields(:,1)==dta(i,2));
           %disp(sprintf('Year %.0f and indxField %.of',dta(i,2),indxField));
           x(indxCategory,1+indxField) = dta(i,3);
       end %for i
       % Add warning to say we should be using the new function
       warning('This is an old function and is superceded by convertNColToMatrix');
    end
    % ======================================================================
    %> @brief converts n col to n-1 dimensional matrix. Last col is the
    %> values and the other cols are the categories
    %> 
    %>
    %> @param dta n col matrix
    %> @param nanVal value to use for Nans [optional]
    %> @retval x matrix, n-1 dims
    %> @retval fields structure with fieldnames field1, field2 etc with a
    %> list in each of the the values representing the corresponding index of
    %> the x matrix
    % ======================================================================
    function [x,fields] = convertNColToMatrix(dta,nanVal)
        if nargin==1
            nanVal = nan;
        end %if
       cats = dta(:,1:end-1);
       catNo = size(dta,2)-1;
       rowNo = size(dta,1);
       fields = struct('field1',sort(unique(dta(:,1))));
       for i=2:catNo
           fname = [sprintf('field%i',i)];
           fields.(fname) = sort(unique(dta(:,i)));
       end %for i
       fieldnm = fieldnames(fields);
       % declare size of x matrix
       switch catNo
           case 2
               x = zeros(size(fields.field1,1),size(fields.field2,1));
           case 3
               x = zeros(size(fields.field1,1),size(fields.field2,1),size(fields.field3,1));
           case 4
               x = zeros(size(fields.field1,1),size(fields.field2,1),...
                   size(fields.field3,1),size(fields.field4,1));
           case 5
               x = zeros(size(fields.field1,1),size(fields.field2,1),...
                   size(fields.field3,1),size(fields.field4,1),size(fields.field5,1));
       end %switch
       for i=1:rowNo
           indxField =struct;
           for j=1:size(fieldnm)
                indxField.(char(fieldnm(j))) = ...
                    find(fields.(char(fieldnm(j)))==dta(i,j));
           end %for j
           % if nan then replace it with nanVal
           if isnan(dta(i,end))
               dta(i,end) = nanVal;
           end %if
           % NOw fill in value
           switch catNo
               case 2
                   x(indxField.field1,indxField.field2) = ...
                       dta(i,end);
               case 3
                   x(indxField.field1,indxField.field2,indxField.field3) = ...
                       dta(i,end);
               case 4
                   x(indxField.field1,indxField.field2,indxField.field3,...
                       indxField.field4) = ...
                       dta(i,end);
               case 5
                   x(indxField.field1,indxField.field2,indxField.field3,...
                       indxField.field4,indxField.field5) = ...
                       dta(i,end);
           end %switch
       end %for i
       
    end
    % ======================================================================
    %> @brief combosFromColumnVector Creates a two col matrix from a column vector. For example,
    %>you may have an array of countries and want to build combinations of
    %>their route distances so you'd need a square matrix of their
    %>combinations
    %>
    %> @param colVector array of values that we want combinations of
    %> @param stripMatches [optional] set equal to 1 if you want to stirp
    %> matching values (like within country trade) or any other number if
    %> not
    %> @retval comboMatrix two col matrix with all combinations of
    %> colVector
    % ======================================================================
    function [comboMatrix] = combosFromColumnVector(colVector,stripMatches)
        if nargin ==1
            stripMatches =1;
        end % if
        importers = repmat(colVector,numel(colVector),1);
        exporters = repmat(colVector',numel(colVector),1);
        exporters = reshape(exporters,numel(exporters),1);
        comboMatrix = [importers exporters];
        if stripMatches ==1
            %strip off rows where first col is equal to second col
            comboMatrix(comboMatrix(:,1)==comboMatrix(:,2),:)=[];
        end %if
    end %function 
    
     % ======================================================================
    %> @brief return parsed string for long and lat from a string to a
    %> numeric value for hte coordinate
    %>
    %> @param string containing lat or lon in format 26-25N or 34-32E
    %> @retval lat or long
    % ======================================================================
    function [coord] = getShellLatLon(LatLon)
        %get array of numeric values
        coord = cell2mat(regexp(LatLon,{'\d'}));
        coord = str2num(LatLon(coord(1:end-2))) + ...
            str2num(LatLon(coord(end-1:end)))/60;
        %if LatLon contains an or a W then multiply by -1
        if ~isempty(regexp(LatLon,'[NW]', 'once'))
            coord = coord * -1;
        end
    end %function getShellLatLon
    %==================================================================
    %> @brief converts an index to a character representing the column in
    %> an excel file
    %>
    %>
    %> @param charIndx
    function [xlsChar]= getXlsChar(charIndx)
        xlsChar = '';
        switch floor(charIndx/26)
            case 0
                xlsChar = '';
            case 1
                xlsChar = 'A';
            case 2

                xlsChar = 'B';
            case 3
                xlsChar = 'C';
            case 4
                xlsChar = 'D';
            case 5
                xlsChar ='E';
        end%switch
        xlsChar = [xlsChar char(96+rem(charIndx,26))];
    end %getXlsChar
   
    %**************************************************
    %> @brief Reads properties from a config file
    %>
    %> @param parameter the parameter to get 
    %> @param fileLoc location of config file
    %>
    %> @retval The parameter value
    function [val]= getConfigProperty(parameter,fileLoc)
        if nargin==1
            fileLoc = 'config.properties';
        end %if
        fileID = fopen(fileLoc);
        data = textscan(fileID, '%s');
        fclose(fileID);
        
        %Find the line we want
        indx = cellfun(@(x) ~isempty(x),regexp(data{1,1},parameter));

        val = data{1,1};
        val = val{indx};
        val = val(regexp(val,'=')+1:end);
        
        
    end %function getConfigProperty
   end % methods Static = true

end
