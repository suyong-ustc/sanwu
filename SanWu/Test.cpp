#include "Test.h"
#include "DIC/DICAlgorithm.h"

using namespace arma;



void SetDICParameters(DICParameters& dic_parameters)
{
	// 感兴趣区域
	dic_parameters.set_roi(50, 450, 50, 450);

	// 网格间距
	dic_parameters.set_grid_step(1);

	// 子区尺寸
	dic_parameters.set_subset_size(29);

	// ZNCC阈值
	dic_parameters.set_zncc_threshold(0.8);

	// 最大迭代次数
	dic_parameters.set_max_iteration_times(10);

	// 误差阈值
	dic_parameters.set_error_threshold(2e-4);

	// 插值阶数
	dic_parameters.set_bspline_interpolation_order(3);

	// 形函数阶数
	dic_parameters.set_shape_function_order(1);
}



bool ReadImage(const std::string& image_path, mat& image)
{
	//std::cout << "Import image with path " << image_path << std::endl;

	// 读取图像
	cv::Mat cvmat = cv::imread(image_path, cv::IMREAD_GRAYSCALE);

	if (cvmat.empty())
	{
		std::cerr << "I can not import the image!" << std::endl;
		return false;
	}

	// 将 opencv 的数据格式转化为 armadillo 矩阵
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
	/***********        读取图像             ************/
	/***************************************************/
	
	const std::string prefix = std::string("..\\images\\Sinusoidal\\T") + std::to_string(period) + "_";
	const std::string refer_image_path = prefix + "0.bmp";
	const std::string deform_image_path = prefix + "1.bmp";

	mat refer_image;
	mat deform_image;
	ReadImage(refer_image_path, refer_image);
	ReadImage(deform_image_path, deform_image);


	/***************************************************/
	/***********        计算参数             ************/
	/***************************************************/

	DICParameters dic_parameters;
	SetDICParameters(dic_parameters);

	// 设置子区大小
	dic_parameters.set_subset_size(subset_size);

	// 设置形函数阶数
	dic_parameters.set_shape_function_order(shape_function_order);

	// 设置权重函数
	mat w = zeros(subset_size, subset_size);

	for (int r = 0; r < subset_size; ++r)
		for (int c = 0; c < subset_size - empty; ++c)
			w(r, c) = 1;

	dic_parameters.set_weight(vectorise(w));


	/***************************************************/
	/***********        相关计算             ************/
	/***************************************************/

	// 网格点
	mat grid_x;
	mat grid_y;
	dic_parameters.grid(grid_x, grid_y);

	// 初始化计算结果
	DICOutput* dic_output = new DICOutput(grid_x, grid_y, ParameterTotalNumber(dic_parameters.shape_function_order()));

	// 设置迭代初值
	dic_output->set_u(sin(2.0 * datum::pi * grid_x / period));

	// 相关计算
	RegisterSubpixelDisplacement(refer_image, deform_image, dic_parameters, dic_output);


	/***************************************************/
	/***********        存储结果             ************/
	/***************************************************/

	// 输出路径
	std::string output_prefix = 
		std::string("..\\results\\Sinusoidal\\T") + std::to_string(period) + 
		"M" + std::to_string(dic_parameters.subset_size()) + 
		"U" + std::to_string(empty) + 
		"N" + std::to_string(dic_parameters.shape_function_order());
	dic_output->write(output_prefix);

	// 释放内存
	delete dic_output;

	std::cout << std::endl;


}





void AnalyzeTransferFunctions()
{
	// 0: 分析不同形函数影响；
	// 1：分析不同子区大小影响；
	// 2：分析不同周期影响
	int mode = 2;

	if (mode == 0)
		std::cout << "Analyze sinusoidal displacement fields corresponding to different shape function orders ..." << std::endl;
	else if (mode == 1)
		std::cout << "Analyze sinusoidal displacement fields corresponding to different subset sizes ..." << std::endl;
	else if (mode == 2)
		std::cout << "Analyze sinusoidal displacement fields corresponding to different periods ..." << std::endl;


	// 参数矩阵
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


	// 相关计算
	for (int i = 0; i < para.n_rows; ++i)
	{
		// 参数
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



void AnalyzeBoundary()
{
	/********************************************/
	/*********         设置参数            *******/
	/********************************************/

	// 初始化计算参数
	DICParameters dic_parameters;
	SetDICParameters(dic_parameters);

	// 设置子区大小
	const int subset_size = 29;
	dic_parameters.set_subset_size(subset_size);

	// 设置形函数阶数
	const int shape_function_order = 2;
	dic_parameters.set_shape_function_order(shape_function_order);

	// 设置感兴趣区域
	const int half_subset_size = (subset_size - 1) / 2;
	const int x0 = 150;
	dic_parameters.set_roi(x0 - half_subset_size, x0 + half_subset_size, 50, 450);


	/********************************************/
	/*********         读取图像            *******/
	/********************************************/

	// 读取图像
	const std::string prefix = std::string("..\\images\\Polynomial\\");
	const std::string refer_image_path = prefix + "a-4n3_0.bmp";
	const std::string deform_image_path = prefix + "a-4n3_1.bmp";

	mat refer_image;
	mat deform_image;
	ReadImage(refer_image_path, refer_image);
	ReadImage(deform_image_path, deform_image);


	/********************************************/
	/*********         相关计算            *******/
	/********************************************/

	// 网格点
	mat grid_x;
	mat grid_y;
	dic_parameters.grid(grid_x, grid_y);

	// 初始化计算结果
	DICOutput* dic_output = new DICOutput(grid_x, grid_y, ParameterTotalNumber(dic_parameters.shape_function_order()));

	// 设定迭代初值
	const mat u = 0.0001 * pow(grid_x - 150, 3);
	const mat v = zeros(grid_x.n_rows, grid_x.n_cols);

	dic_output->set_u(u);
	dic_output->set_v(v);

	// 相关计算
	RegisterSubpixelDisplacementWithBoundaryCheck(refer_image, deform_image, dic_parameters, dic_output);

	// 输出结果
	//std::string output_prefix = std::string("..\\results\\Sinusoidal\\T") + std::to_string(period)
	//	+ "N" + std::to_string(dic_parameters.shape_function_order())
	//	+ "D" + std::to_string(nul)
	//	+ "M" + std::to_string(dic_parameters.subset_size()) + "_";
	dic_output->write("n3");

	// 释放内存
	delete dic_output;

	std::cout << "\nJob Done!" << std::endl;

}
