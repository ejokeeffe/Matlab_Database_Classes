%> @file DataBase.m
%> @brief Abstract class for inheritance for DataBase subclasses (eg for
%xls, access, mysql etc)

%> @mainpage Overview

%> @section intro Introduction

%> Authors: Eoin O'Keeffe,

%> Date initiated: 13/8/10

%> Version: 1.0

%> To ensure consistency and extensibility, this class is inherited by all
%database classes
classdef DataBaseAbs < handle
    %DATABASE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        db = 'inputs\trade\Routes.xls';
        root = '';
    end
    methods (Sealed = true)
    end % methods
    methods (Abstract = true)
         success = open(obj)
        success = close(obj)
         x = rows(obj,tableName) %Returns number of rows
         x    = columns(obj,tableName) %Returns number of columns
         x = getAll(obj,tableName) %array (for excel a cell-array is passed)
         x = getCol(obj,tableName,index) %array (for excel a cell-array is passed)
         x = getRow(obj,tableName,index) %array (for excel a cell-array is passed)
         x = getOne(obj,tableName,index) %returns a single value
         success = setAll(obj,tableName,array) 
         success = setCol(obj,tableName,x,array)
         success = setRow(obj,tableName,y,array)
         success = setOne(obj,tableName,x,y, value)
    end
    
end

