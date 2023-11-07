#include <iostream>
#include <QStringList>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <armadillo>
#include "DIC/DICAlgorithm.h"
#include "DIC/DICParameters.h"
using namespace arma;



void SetDICParameters(DICParameters& dic_parameters)
{
	// ����Ȥ����
	dic_parameters.set_roi(50, 450, 50, 450);
	
	// ������
	dic_parameters.set_grid_step(1);

	// �����ߴ�
	dic_parameters.set_subset_size(29);

	// ZNCC��ֵ
	dic_parameters.set_zncc_threshold(0.8);

	// ����������
	dic_parameters.set_max_iteration_times(10);

	// �����ֵ
	dic_parameters.set_error_threshold(2e-4);

	// ��ֵ����
	dic_parameters.set_bspline_interpolation_order(3);

	// �κ�������
	dic_parameters.set_shape_function_order(1);
}



bool ReadImage(const std::string& image_path, mat& image)
{
	std::cout << "Import image with path " << image_path << std::endl;

	// ��ȡͼ��
	cv::Mat cvmat = cv::imread(image_path, cv::IMREAD_GRAYSCALE);

	if (cvmat.empty())
	{
		std::cerr << "I can not import the image!" << std::endl;
		return false;
	}

	// �� opencv �����ݸ�ʽת��Ϊ armadillo ����
	image.zeros(cvmat.rows, cvmat.cols);

	for (int r = 0; r < cvmat.rows; ++r)
	{
		for (int c = 0; c < cvmat.cols; ++c)
		{
			image(r, c) = cvmat.at<uchar>(r, c);
		}

	}

	return true;
}



void AnalyzeSinusoidalDisplacement()
{
	// 0: ������ͬ�κ���Ӱ�죻
	// 1��������ͬ������СӰ�죻
	// 2��������ͬ����Ӱ��
	int mode = 1;

	if (mode == 0)
		std::cout << "Analyze sinusoidal displacement fields corresponding to different shape function orders ..." << std::endl;
	else if (mode == 1)
		std::cout << "Analyze sinusoidal displacement fields corresponding to different subset sizes ..." << std::endl;
	else if (mode == 2)
		std::cout << "Analyze sinusoidal displacement fields corresponding to different periods ..." << std::endl;


	// ��������
	imat para;

	if (mode == 0)
	{
		const int T = 20;
		const int subset_size = 15;
		const ivec shape_function_order = regspace<ivec>(0, 1, 5);

		para.zeros(shape_function_order.n_elem, 3);
		for (int i = 0; i < shape_function_order.n_elem; ++i)
		{
			para(i, 0) = T;
			para(i, 1) = subset_size;
			para(i, 2) = shape_function_order(i);
		}
	}
	else if (mode == 1)
	{
		const int T = 40;
		const ivec subset_size = regspace<ivec>(19, -2, 9);
		const int shape_function_order = 4;

		para.zeros(subset_size.n_elem, 3);
		for (int i = 0; i < subset_size.n_elem; ++i)
		{
			para(i, 0) = T;
			para(i, 1) = subset_size(i);
			para(i, 2) = shape_function_order;
		}
	}
	else if (mode == 2)
	{
		const ivec T = regspace<ivec>(8, -4, 4);
		const int subset_size = 29;
		const int shape_function_order = 4;

		para.zeros(T.n_elem, 3);
		for (int i = 0; i < T.n_elem; ++i)
		{
			para(i, 0) = T(i);
			para(i, 1) = subset_size;
			para(i, 2) = shape_function_order;
		}

	}


	// ��ʼ���������
	DICParameters dic_parameters;
	SetDICParameters(dic_parameters);

	// ��ؼ���
	for (int i = 0; i < para.n_rows; ++i)
	{
		// ����
		const int period = para(i, 0);
		const int subset_size = para(i, 1);
		const int shape_function_order = para(i, 2);

		std::cout << "Start analyze period: " << period << ";\tsubset size: " << subset_size << ";\tshape function order: " << shape_function_order << std::endl;

		// ��ȡͼ��
		const std::string prefix = std::string("..\\images\\Sinusoidal\\T") + std::to_string(period) + "_";
		const std::string refer_image_path = prefix + "0.bmp";
		const std::string deform_image_path = prefix + "1.bmp";

		mat refer_image;
		mat deform_image;
		ReadImage(refer_image_path, refer_image);
		ReadImage(deform_image_path, deform_image);

		// ����������С
		dic_parameters.set_subset_size(subset_size);


		// �����κ�������
		dic_parameters.set_shape_function_order(shape_function_order);

		// �����
		mat grid_x;
		mat grid_y;
		dic_parameters.grid(grid_x, grid_y);

		// ��ʼ��������
		DICOutput* dic_output = new DICOutput(grid_x, grid_y, ParameterTotalNumber(dic_parameters.shape_function_order()));

		// ���Ƴ�ֵ
		EstimateInitialDisplacement(refer_image, deform_image, dic_parameters, dic_output);

		if (dic_parameters.subset_size() < 20 && dic_parameters.shape_function_order() == 4)
		{
			mat u = sin(2 * datum::pi * grid_x / period);
			dic_output->set_u(u);
		}


		// ��ؼ���
		RegisterSubpixelDisplacement(refer_image, deform_image, dic_parameters, dic_output);

		// ������
		std::string output_prefix = std::string("..\\results\\Sinusoidal\\T") + std::to_string(period) + "N" + std::to_string(dic_parameters.shape_function_order()) + "M" + std::to_string(dic_parameters.subset_size()) + "_";
		dic_output->write(output_prefix);

		// �ͷ��ڴ�
		delete dic_output;

		std::cout << std::endl;

	}

	std::cout << "Job Done!" << std::endl;
}



void AnalyzeGaussianDisplacement()
{
	// 0: ������ͬ�κ���Ӱ�죻
	// 1��������ͬ������СӰ�죻
	// 2��������ͬ����Ӱ��
	int mode = 1;

	if (mode == 0)
		std::cout << "Analyze gaussian displacement fields corresponding to different shape function orders ..." << std::endl;
	else if (mode == 1)
		std::cout << "Analyze gaussian displacement fields corresponding to different subset sizes ..." << std::endl;
	else if (mode == 2)
		std::cout << "Analyze gaussian displacement fields corresponding to different widths ..." << std::endl;


	// ��������
	imat para;

	if (mode == 0)
	{
		const int C = 10;
		const int subset_size = 29;
		const ivec shape_function_order = regspace<ivec>(0, 1, 5);

		para.zeros(shape_function_order.n_elem, 3);
		for (int i = 0; i < shape_function_order.n_elem; ++i)
		{
			para(i, 0) = C;
			para(i, 1) = subset_size;
			para(i, 2) = shape_function_order(i);
		}
	}
	else if (mode == 1)
	{
		const int C = 10;
		const ivec subset_size = regspace<ivec>(19, -2, 9);
		const int shape_function_order = 4;

		para.zeros(subset_size.n_elem, 3);
		for (int i = 0; i < subset_size.n_elem; ++i)
		{
			para(i, 0) = C;
			para(i, 1) = subset_size(i);
			para(i, 2) = shape_function_order;
		}
	}
	else if (mode == 2)
	{
		const ivec C = regspace<ivec>(50, -2, 2);
		const int subset_size = 29;
		const int shape_function_order = 4;

		para.zeros(C.n_elem, 3);
		for (int i = 0; i < C.n_elem; ++i)
		{
			para(i, 0) = C(i);
			para(i, 1) = subset_size;
			para(i, 2) = shape_function_order;
		}

	}


	// ��ʼ���������
	DICParameters dic_parameters;
	SetDICParameters(dic_parameters);

	// ��ؼ���
	for (int i = 0; i < para.n_rows; ++i)
	{
		// ����
		const int width = para(i, 0);
		const int subset_size = para(i, 1);
		const int shape_function_order = para(i, 2);

		std::cout << "Start analyze width: " << width << ";\tsubset size: " << subset_size << ";\tshape function order: " << shape_function_order << std::endl;

		// ��ȡͼ��
		const std::string prefix = std::string("..\\images\\Gaussian\\C") + std::to_string(width) + "_";
		const std::string refer_image_path = prefix + "0.bmp";
		const std::string deform_image_path = prefix + "1.bmp";

		mat refer_image;
		mat deform_image;
		ReadImage(refer_image_path, refer_image);
		ReadImage(deform_image_path, deform_image);

		// ����������С
		dic_parameters.set_subset_size(subset_size);

		// �����κ�������
		dic_parameters.set_shape_function_order(shape_function_order);

		// �����
		mat grid_x;
		mat grid_y;
		dic_parameters.grid(grid_x, grid_y);

		// ��ʼ��������
		DICOutput* dic_output = new DICOutput(grid_x, grid_y, ParameterTotalNumber(dic_parameters.shape_function_order()));

		// ���Ƴ�ֵ
		EstimateInitialDisplacement(refer_image, deform_image, dic_parameters, dic_output);

		if (dic_parameters.subset_size() < 20 && dic_parameters.shape_function_order() == 4)
		{			
			mat t = (dic_output->x() - 250) / width;
			mat u = exp(-t % t);
			dic_output->set_u(u);
		}

		// ��ؼ���
		RegisterSubpixelDisplacement(refer_image, deform_image, dic_parameters, dic_output);

		// ������
		std::string output_prefix = std::string("..\\results\\Gaussian\\C") + std::to_string(width) + "N" + std::to_string(dic_parameters.shape_function_order()) + "M" + std::to_string(dic_parameters.subset_size()) + "_";
		dic_output->write(output_prefix);

		// �ͷ��ڴ�
		delete dic_output;

		std::cout << std::endl;

	}

	std::cout << "Job Done!" << std::endl;
}



void AnalyzeTransferFunctionOfUnfullSubset()
{
	// ����
	const int period = 40;
	const int subset_size = 19;
	const int shape_function_order = 1;

	// Ȩ��
	const int nul = 3;

	mat w = zeros(subset_size, subset_size);
	for (int r = 0; r < subset_size; ++r)
		for (int c = 0; c < subset_size - nul; ++c)
			w(r, c) = 1;

	const vec weight = vectorise(w);

	std::cout << "Start analyze period: " << period
		<< ";\tsubset size: " << subset_size
		<< ";\tnul: " << nul
		<< ";\tshape function order: " << shape_function_order << std::endl;

	// ��ȡͼ��
	const std::string prefix = std::string("..\\images\\Sinusoidal\\T") + std::to_string(period) + "_";
	const std::string refer_image_path = prefix + "0.bmp";
	const std::string deform_image_path = prefix + "1.bmp";

	mat refer_image;
	mat deform_image;
	ReadImage(refer_image_path, refer_image);
	ReadImage(deform_image_path, deform_image);

	// ��ʼ���������
	DICParameters dic_parameters;
	SetDICParameters(dic_parameters);

	// ����������С
	dic_parameters.set_subset_size(subset_size);

	// �����κ�������
	dic_parameters.set_shape_function_order(shape_function_order);

	// ����Ȩ�غ���
	dic_parameters.set_weight(weight);

	// �����
	mat grid_x;
	mat grid_y;
	dic_parameters.grid(grid_x, grid_y);

	// ��ʼ��������
	DICOutput* dic_output = new DICOutput(grid_x, grid_y, ParameterTotalNumber(dic_parameters.shape_function_order()));

	// ���Ƴ�ֵ
	EstimateInitialDisplacement(refer_image, deform_image, dic_parameters, dic_output);

	if (dic_parameters.subset_size() < 20 && dic_parameters.shape_function_order() == 4)
	{
		mat u = sin(2 * datum::pi * grid_x / period);
		dic_output->set_u(u);
	}


	// ��ؼ���
	RegisterSubpixelDisplacement(refer_image, deform_image, dic_parameters, dic_output);

	// ������
	std::string output_prefix = std::string("..\\results\\Sinusoidal\\T") + std::to_string(period)
		+ "N" + std::to_string(dic_parameters.shape_function_order())
		+ "D" + std::to_string(nul)
		+ "M" + std::to_string(dic_parameters.subset_size()) + "_";
	dic_output->write(output_prefix);

	// �ͷ��ڴ�
	delete dic_output;

	std::cout << "\nJob Done!" << std::endl;

}



void AnalyzeBoundary()
{
	/********************************************/
	/*********         ���ò���            *******/
	/********************************************/

	// ��ʼ���������
	DICParameters dic_parameters;
	SetDICParameters(dic_parameters);

	// ����������С
	const int subset_size = 19;
	dic_parameters.set_subset_size(subset_size);

	// �����κ�������
	const int shape_function_order = 1;
	dic_parameters.set_shape_function_order(shape_function_order);

	// ���ø���Ȥ����
	const int half_subset_size = (subset_size - 1) / 2;
	const int x0 = 250;
	dic_parameters.set_roi(x0 - 20, x0 + 20, 50, 450);


	/********************************************/
	/*********         ��ȡͼ��            *******/
	/********************************************/

	// ��ȡͼ��
	const std::string prefix = std::string("..\\images\\Polynomial\\");
	const std::string refer_image_path = prefix + "a-3n3_0.bmp";
	const std::string deform_image_path = prefix + "a-3n3_1.bmp";

	mat refer_image;
	mat deform_image;
	ReadImage(refer_image_path, refer_image);
	ReadImage(deform_image_path, deform_image);


	/********************************************/
	/*********         ��ȡͼ��            *******/
	/********************************************/

	// �����
	mat grid_x;
	mat grid_y;
	dic_parameters.grid(grid_x, grid_y);

	// ��ʼ��������
	DICOutput* dic_output = new DICOutput(grid_x, grid_y, ParameterTotalNumber(dic_parameters.shape_function_order()));

	// ���Ƴ�ֵ
	EstimateInitialDisplacement(refer_image, deform_image, dic_parameters, dic_output);

	// ��ؼ���
	RegisterSubpixelDisplacement(refer_image, deform_image, dic_parameters, dic_output);

	// ������
	//std::string output_prefix = std::string("..\\results\\Sinusoidal\\T") + std::to_string(period)
	//	+ "N" + std::to_string(dic_parameters.shape_function_order())
	//	+ "D" + std::to_string(nul)
	//	+ "M" + std::to_string(dic_parameters.subset_size()) + "_";
	dic_output->write("n3");

	// �ͷ��ڴ�
	delete dic_output;

	std::cout << "\nJob Done!" << std::endl;

}



int main()
{
	AnalyzeBoundary();
}