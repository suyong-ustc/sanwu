#include "Test.h"
#include "DIC/DICAlgorithm.h"

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
	//std::cout << "Import image with path " << image_path << std::endl;

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


bool TestNormalizeVectorize()
{
	mat a = randu(3, 3);

	mat w = zeros(3, 3);
	for (int r = 1; r < 3; ++r)
		for (int c = 1; c < 3; ++c)
			w(r, c) = 1;

	const double zncc = ZNCC(a, a, w);
	std::cout << "the correaltion function is " << zncc << std::endl;

	return true;
}



void AnalyzeTransferFunction(const int& period, const int& subset_size, const int& empty, const int& shape_function_order)
{
	std::cout << "Start analyze period: " << period << ";\t"
		<< "subset size: " << subset_size << ";\t"
		<< "empty: " << empty << ";\t"
		<< "shape function order: " << shape_function_order << std::endl;

	/***************************************************/
	/***********        ��ȡͼ��             ************/
	/***************************************************/
	
	const std::string prefix = std::string("..\\images\\Sinusoidal\\T") + std::to_string(period) + "_";
	const std::string refer_image_path = prefix + "0.bmp";
	const std::string deform_image_path = prefix + "1.bmp";

	mat refer_image;
	mat deform_image;
	ReadImage(refer_image_path, refer_image);
	ReadImage(deform_image_path, deform_image);


	/***************************************************/
	/***********        �������             ************/
	/***************************************************/

	DICParameters dic_parameters;
	SetDICParameters(dic_parameters);

	// ����������С
	dic_parameters.set_subset_size(subset_size);

	// �����κ�������
	dic_parameters.set_shape_function_order(shape_function_order);

	// ����Ȩ�غ���
	mat w = zeros(subset_size, subset_size);

	for (int r = 0; r < subset_size; ++r)
		for (int c = 0; c < subset_size - empty; ++c)
			w(r, c) = 1;

	dic_parameters.set_weight(vectorise(w));


	/***************************************************/
	/***********        ��ؼ���             ************/
	/***************************************************/

	// �����
	mat grid_x;
	mat grid_y;
	dic_parameters.grid(grid_x, grid_y);

	// ��ʼ��������
	DICOutput* dic_output = new DICOutput(grid_x, grid_y, ParameterTotalNumber(dic_parameters.shape_function_order()));

	// ���õ�����ֵ
	dic_output->set_u(sin(2.0 * datum::pi * grid_x / period));

	// ��ؼ���
	RegisterSubpixelDisplacement(refer_image, deform_image, dic_parameters, dic_output);


	/***************************************************/
	/***********        �洢���             ************/
	/***************************************************/

	// ���·��
	std::string output_prefix = 
		std::string("..\\results\\Sinusoidal\\T") + std::to_string(period) + 
		"M" + std::to_string(dic_parameters.subset_size()) + 
		"U" + std::to_string(empty) + 
		"N" + std::to_string(dic_parameters.shape_function_order());
	dic_output->write(output_prefix);

	// �ͷ��ڴ�
	delete dic_output;

	std::cout << std::endl;


}





void AnalyzeTransferFunctions()
{
	// 0: ������ͬ�κ���Ӱ�죻
	// 1��������ͬ������СӰ�죻
	// 2��������ͬ����Ӱ��
	int mode = 2;

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
		const int period = 50;
		const int subset_size = 29;
		const ivec shape_function_order = regspace<ivec>(0, 1, 2);

		para.zeros(shape_function_order.n_elem, 3);
		for (int i = 0; i < shape_function_order.n_elem; ++i)
		{
			para(i, 0) = period;
			para(i, 1) = subset_size;
			para(i, 2) = shape_function_order(i);
		}
	}
	else if (mode == 1)
	{
		const int period = 50;
		const ivec subset_size = regspace<ivec>(19, -10, 19);
		const int shape_function_order = 2;

		para.zeros(subset_size.n_elem, 3);
		for (int i = 0; i < subset_size.n_elem; ++i)
		{
			para(i, 0) = period;
			para(i, 1) = subset_size(i);
			para(i, 2) = shape_function_order;
		}
	}
	else if (mode == 2)
	{
		const ivec period = regspace<ivec>(60, -10, 60);
		const int subset_size = 29;
		const int shape_function_order = 2;

		para.zeros(period.n_elem, 3);
		for (int i = 0; i < period.n_elem; ++i)
		{
			para(i, 0) = period(i);
			para(i, 1) = subset_size;
			para(i, 2) = shape_function_order;
		}

	}


	// ��ؼ���
	for (int i = 0; i < para.n_rows; ++i)
	{
		// ����
		const int period = para(i, 0);
		const int subset_size = para(i, 1);
		const int shape_function_order = para(i, 2);

		for (int empty = 0; empty < 0.75 * subset_size; ++empty)
		{
			AnalyzeTransferFunction(period, subset_size, empty, shape_function_order);
		}

	}

	std::cout << "Job Done!" << std::endl;

}



void AnalyzeBoundary(const int& polynomial_order, const int& subset_size, const int& shape_function_order)
{
	std::cout << "Start analyze polynomial order: " << polynomial_order << ";\t"
		<< "subset size: " << subset_size << ";\t"
		<< "shape function order: " << shape_function_order << std::endl;

	/********************************************/
	/*********         ���ò���            *******/
	/********************************************/

	// ��ʼ���������
	DICParameters dic_parameters;
	SetDICParameters(dic_parameters);

	// ����������С
	dic_parameters.set_subset_size(subset_size);

	// �����κ�������
	dic_parameters.set_shape_function_order(shape_function_order);

	// ���ø���Ȥ����
	if (polynomial_order != 4) {
		const int half_subset_size = (subset_size - 1) / 2;
		const int x0 = 150;
		dic_parameters.set_roi(x0 - half_subset_size, x0 + half_subset_size, 50, 450);
	}
	else if (polynomial_order == 4) {
		const int half_subset_size = (subset_size - 1) / 2;
		const int x0 = 40;
		dic_parameters.set_roi(x0 - half_subset_size, x0 + half_subset_size, 50, 450);
	}


	/********************************************/
	/*********         ��ȡͼ��            *******/
	/********************************************/

	// ��ȡͼ��
	const std::string prefix = std::string("..\\images\\Polynomial\\");
	std::string refer_image_path;
	std::string deform_image_path;
	if (polynomial_order == 2) {
		refer_image_path = prefix + "A-3L2_0.bmp";
		deform_image_path = prefix + "A-3L2_1.bmp";
	}
	else if (polynomial_order == 3) {
		refer_image_path = prefix + "A-4L3_0.bmp";
		deform_image_path = prefix + "A-4L3_1.bmp";
	}
	else if (polynomial_order == 4) {
		refer_image_path = prefix + "A-6L4_0.bmp";
		deform_image_path = prefix + "A-6L4_1.bmp";
	}
	else if (polynomial_order == 5) {
		refer_image_path = prefix + "A-7L5_0.bmp";
		deform_image_path = prefix + "A-7L5_1.bmp";
	}

	mat refer_image;
	mat deform_image;
	ReadImage(refer_image_path, refer_image);
	ReadImage(deform_image_path, deform_image);


	/********************************************/
	/*********         ��ؼ���            *******/
	/********************************************/

	// �����
	mat grid_x;
	mat grid_y;
	dic_parameters.grid(grid_x, grid_y);

	// ��ʼ��������
	DICOutput* dic_output = new DICOutput(grid_x, grid_y, ParameterTotalNumber(dic_parameters.shape_function_order()));

	// �趨������ֵ
	const mat v = zeros(grid_x.n_rows, grid_x.n_cols);
	mat u;
	if (polynomial_order == 2)
		u = 1e-3 * pow(grid_x - 150, 2);
	else if (polynomial_order == 3)
		u = 1e-4 * pow(grid_x - 150, 3);
	else if (polynomial_order == 4)
		u = 1e-6 * pow(grid_x - 40, 4);
	else if (polynomial_order == 5)
		u = 1e-7 * pow(grid_x - 150, 5);

	dic_output->set_u(u);
	dic_output->set_v(v);

	// ��ؼ���
	RegisterSubpixelDisplacementWithBoundaryCheck(refer_image, deform_image, dic_parameters, dic_output);

	// ������
	std::string output_prefix =
		std::string("..\\results\\Polynomial\\L") + std::to_string(polynomial_order)
		+ "N" + std::to_string(shape_function_order)
		+ "M" + std::to_string(subset_size);

	dic_output->write(output_prefix);

	// �ͷ��ڴ�
	delete dic_output;

	std::cout << "\nJob Done!" << std::endl;

}


void AnalyzeBoundaries()
{
	// 0: ������ͬ����ʽ����Ӱ�죻
	// 1��������ͬ������СӰ�죻

	int mode = 1;

	if (mode == 0)
		std::cout << "Analyze boundary displacement fields corresponding to different polynomial orders ..." << std::endl;
	else if (mode == 1)
		std::cout << "Analyze sinusoidal displacement fields corresponding to different subset sizes ..." << std::endl;


	// ��������
	imat para;

	if (mode == 0)
	{
		const ivec polynomail_order = regspace<ivec>(2, 1, 5);
		const int subset_size = 29;
		const int shape_function_order = 2;

		para.zeros(polynomail_order.n_elem, 3);
		for (int i = 0; i < polynomail_order.n_elem; ++i)
		{
			para(i, 0) = polynomail_order(i);
			para(i, 1) = subset_size;
			para(i, 2) = shape_function_order;
		}
	}
	else if (mode == 1)
	{
		const int polynomail_order = 3;
		const ivec subset_size = regspace<ivec>(9, 10, 9);
		const int shape_function_order = 0;

		para.zeros(subset_size.n_elem, 3);
		for (int i = 0; i < subset_size.n_elem; ++i)
		{
			para(i, 0) = polynomail_order;
			para(i, 1) = subset_size(i);
			para(i, 2) = shape_function_order;
		}
	}


	// ��ؼ���
	for (int i = 0; i < para.n_rows; ++i)
	{
		// ����
		const int polynomial_order = para(i, 0);
		const int subset_size = para(i, 1);
		const int shape_function_order = para(i, 2);

		AnalyzeBoundary(polynomial_order, subset_size, shape_function_order);

	}

	std::cout << "Job Finished!" << std::endl;


}
