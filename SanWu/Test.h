#pragma once
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <armadillo>
#include "DIC/DICParameters.h"



void SetDICParameters(DICParameters& dic_parameters);



bool ReadImage(const std::string& image_path, arma::mat& image);



bool TestNormalizeVectorize();



void AnalyzeTransferFunction(const int& period, const int& subset_size, const int& empty, const int& shape_function_order);



void AnalyzeTransferFunctions();





void AnalyzeBoundary();



void AnalyzeBoundaries();