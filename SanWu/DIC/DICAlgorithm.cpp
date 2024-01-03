#include "DICAlgorithm.h"
#include "..\Interpolation\BicubicMOMSInterpolator.h"
#include "..\Interpolation\BiquinticBSplineInterpolatror.h"
#include "..\Interpolation\BisepticBSplineInterpolator.h"
using namespace arma;




/********************************************************************************************************/
/******************************                    形函数              **********************************/
/********************************************************************************************************/


bool ShapeFunction(const double& x0, const double& y0, const mat& x, const mat& y, const vec& p, const int& order, mat& deform_x, mat& deform_y)
{
	bool ok = false;

	if (order == 0)
		ok = ZeroOrderShapeFunction(x0, y0, x, y, p, deform_x, deform_y);
	else if (order == 1)
		ok = FirstOrderShapeFunction(x0, y0, x, y, p, deform_x, deform_y);
	else if (order == 2)
		ok = SecondOrderShapeFunction(x0, y0, x, y, p, deform_x, deform_y);

	return ok;
}



bool ZeroOrderShapeFunction(const double& x0, const double& y0, const mat& x, const mat& y, const vec& p, mat& deform_x, mat& deform_y)
{
	// 零阶形函数参数
	const double u = p(0);
	const double v = p(1);

	// 根据零阶形函数估计变形后子区位置
	deform_x = x + u + x0;
	deform_y = y + v + y0;

	return true;
}



bool FirstOrderShapeFunction(const double& x0, const double& y0, const mat& x, const mat& y, const vec& p, mat& deform_x, mat& deform_y)
{
	// 零阶形函数参数
	const double u = p(0);
	const double v = p(1);

	// 根据零阶形函数估计变形后子区位置
	deform_x = x + u + x0;
	deform_y = y + v + y0;

	// 一阶形函数参数
	const double ux = p(2);
	const double uy = p(3);

	const double vx = p(4);
	const double vy = p(5);

	// 根据一阶形函数估计变形后子区位置
	deform_x = deform_x + ux * x + uy * y;
	deform_y = deform_y + vx * x + vy * y;

	return true;
}



bool SecondOrderShapeFunction(const double& x0, const double& y0, const mat& x, const mat& y, const vec& p, mat& deform_x, mat& deform_y)
{
	// 零阶形函数参数
	const double u = p(0);
	const double v = p(1);

	// 根据零阶形函数估计变形后子区位置
	deform_x = x + u + x0;
	deform_y = y + v + y0;

	// 一阶形函数参数
	const double ux = p(2);
	const double uy = p(3);

	const double vx = p(4);
	const double vy = p(5);

	// 根据一阶形函数估计变形后子区位置
	deform_x = deform_x + ux * x + uy * y;
	deform_y = deform_y + vx * x + vy * y;

	// 二阶形函数参数
	const double uxx = p(6);
	const double uxy = p(7);
	const double uyy = p(8);

	const double vxx = p(9);
	const double vxy = p(10);
	const double vyy = p(11);

	// 根据形函数估计变形后子区位置
	const mat xx = x % x;
	const mat xy = x % y;
	const mat yy = y % y;

	deform_x = deform_x + uxx * xx + uxy * xy + uyy * yy;
	deform_y = deform_y + vxx * xx + vxy * xy + vyy * yy;

	return true;
}




/*******************************************************************************************************************/
/**********************************               伪逆矩阵                ********************************************/
/******************************************************************************************************************/



bool PseudoInverseMatrix(const mat& gx, const mat& gy, const mat& x, const mat& y, const vec& w, const int& order, mat& pseudo)
{
	bool ok = false;

	if (order == 0)
		ok = ZeroOrderPseudoInverseMatrix(gx, gy, x, y, w, pseudo);
	else if (order == 1)
		ok = FirstOrderPseudoInverseMatrix(gx, gy, x, y, w, pseudo);
	else if (order == 2)
		ok = SecondOrderPseudoInverseMatrix(gx, gy, x, y, w, pseudo);

	return ok;
}



bool ZeroOrderPseudoInverseMatrix(const mat& gx, const mat& gy, const mat& x, const mat& y, const vec& weight, mat& pseudo)
{
	// 雅可比矩阵
	mat jacobian(gx.n_elem, 2);

	jacobian.col(0) = vectorise(gx);
	jacobian.col(1) = vectorise(gy);

	// 伪逆矩阵
	const mat wj = jacobian.each_col() % weight;

	mat hessian(2, 2);
	for (int r = 0; r < 2; ++r)
		for (int c = 0; c < 2; ++c)
			hessian(r, c) = dot(jacobian.col(r), wj.col(c));

	pseudo = inv_sympd(hessian) * wj.t();

	return true;
}



bool FirstOrderPseudoInverseMatrix(const mat& gx, const mat& gy, const mat& x, const mat& y, const vec& weight, mat& pseudo)
{
	// 雅可比矩阵
	mat jacobian(gx.n_elem, 6);

	jacobian.col(0) = vectorise(gx);
	jacobian.col(1) = vectorise(gy);

	jacobian.col(2) = vectorise(gx % x);
	jacobian.col(3) = vectorise(gx % y);
	jacobian.col(4) = vectorise(gy % x);
	jacobian.col(5) = vectorise(gy % y);

	// 伪逆矩阵
	const mat wj = jacobian.each_col() % weight;

	mat hessian(6, 6);
	for (int r = 0; r < 6; ++r)
		for (int c = 0; c < 6; ++c)
			hessian(r, c) = dot(jacobian.col(r), wj.col(c));

	pseudo = inv_sympd(hessian) * wj.t();

	return true;
}



bool SecondOrderPseudoInverseMatrix(const mat& gx, const mat& gy, const mat& x, const mat& y, const vec& weight, mat& pseudo)
{
	// 雅可比矩阵
	mat jacobian = zeros<mat>(gx.n_elem, 12);

	jacobian.col(0) = vectorise(gx);
	jacobian.col(1) = vectorise(gy);

	jacobian.col(2) = vectorise(gx % x);
	jacobian.col(3) = vectorise(gx % y);
	jacobian.col(4) = vectorise(gy % x);
	jacobian.col(5) = vectorise(gy % y);

	jacobian.col(6) = vectorise(gx % x % x);
	jacobian.col(7) = vectorise(gx % x % y);
	jacobian.col(8) = vectorise(gx % y % y);
	jacobian.col(9) = vectorise(gy % x % x);
	jacobian.col(10) = vectorise(gy % x % y);
	jacobian.col(11) = vectorise(gy % y % y);

	// 伪逆矩阵
	const mat wj = jacobian.each_col() % weight;

	mat hessian(12, 12);
	for (int r = 0; r < 12; ++r)
		for (int c = 0; c < 12; ++c)
			hessian(r, c) = dot(jacobian.col(r), wj.col(c));

	pseudo = inv_sympd(hessian) * wj.t();

	return true;
}





/*******************************************************************************************************************/
/**********************************            计算全场变形             ********************************************/
/******************************************************************************************************************/


DICOutput* RigsterFullFieldDisplacement(const arma::mat& refer_image, const arma::mat& deform_image, const DICParameters& dic_parameters)
{
	// 网格点
	mat grid_x;
	mat grid_y;
	dic_parameters.grid(grid_x, grid_y);

	// 初始化计算结果
	DICOutput* dic_output = new DICOutput(grid_x, grid_y, ParameterTotalNumber(dic_parameters.shape_function_order()));

	// 估计初值
	EstimateInitialDisplacement(refer_image, deform_image, dic_parameters, dic_output);

	// 相关计算
	RegisterSubpixelDisplacement(refer_image, deform_image, dic_parameters, dic_output);


	return dic_output;
}



bool EstimateInitialDisplacement(const mat& refer_image, const mat& deform_image, const DICParameters& dic_parameters, DICOutput* dic_output)
{
	// 网格点
	const mat x = dic_output->x();
	const mat y = dic_output->y();

	// 初始化
	mat u(x.n_rows, x.n_cols, fill::zeros);
	mat v(x.n_rows, x.n_cols, fill::zeros);

	// 赋初值
	//bool is_harmonic = true;

	//if (is_harmonic)
	//{
	//	double a = 1;
	//	double T = 40;
	//	double b = 0;

	//	u = a * sin(2 * datum::pi * x / T);
	//	//u = a * sin(2 * datum::pi * x / T) % sin(2 * datum::pi * y / T);
	//}
	//else
	//{
	//	double a = 1;
	//	double x0 = 250;
	//	double c = 40;

	//	mat t = (x - x0) / c;
	//	u = a * exp(-t % t);
	//}

	//const double a = 0.001;
	//u = a * pow(x - 250, 3);


	dic_output->set_u(u);
	dic_output->set_v(v);

	return true;
}



bool RegisterSubpixelDisplacement(const mat& refer_image, const mat& deform_image, const DICParameters& dic_parameters, DICOutput* dic_output)
{
	// 网格点
	const mat x = dic_output->x();
	const mat y = dic_output->y();

	// "子区点"相对于"子区中心"的坐标
	mat delta_x, delta_y;
	SubsetPointsRelativeCoordinate(dic_parameters.subset_size(), delta_x, delta_y);

	// 变形图插值系数
	Interpolator* deform_image_interpolator = ConstructInterpolator(deform_image, dic_parameters.bspline_interpolation_order());
	
	int num_calculated = 0;
	const int tick = x.n_elem / 10;

#	pragma omp parallel for
	for (int r = 0; r < x.n_rows; ++r)
	{
		for (int c = 0; c < x.n_cols; ++c)
		{
			// 显示计算比例
#			pragma omp critical
			++num_calculated;
			if (num_calculated % tick == 0)
				printf("%.2f%%\t", 100. * num_calculated / x.n_elem);

			// 参考子区灰度
			const mat refer_subset_intensities = SubsetIntensities(refer_image, x(r, c), y(r,c), dic_parameters.half_subset_size());

			// 将参考子区灰度归一化
			vec normalized_refer_subset_intensities;
			NormalizeVectorize(refer_subset_intensities, normalized_refer_subset_intensities);

			// 迭代初值
			vec warp_function(dic_output->warp_function(r, c));

			// 迭代次数
			uword iter = 0;

			mat deform_x, deform_y;			// 变形子区点的图像坐标
			mat deform_subset_intensities;	// 变形子区点的灰度
			mat deform_subset_gradient_x;	// 变形子区点的灰度梯度
			mat deform_subset_gradient_y;
			mat pseudo;						// 伪逆矩阵
			vec normalized_deform_subset_intensities;	// 归一化的变形子区灰度

			// 迭代
			vec dp = ones<vec>(warp_function.n_elem);	// 迭代增量
			while (!isConverge(dp,dic_parameters.error_threshold()) && iter <= dic_parameters.max_iteration_times())
			{
				// 确定变形子区点的"图像坐标"
				ShapeFunction(x(r, c), y(r, c), delta_x, delta_y, warp_function, dic_parameters.shape_function_order(), deform_x, deform_y);

				//// 判断是否越界
				//if (deform_x.min() < 10 || deform_x.max() > 490 || deform_y.min() < 10 || deform_y.max() > 490)
				//	break;

				// 变形子区点的灰度
				deform_subset_intensities = deform_image_interpolator->Values(deform_x, deform_y);
				const double gg = NormalizeVectorize(deform_subset_intensities, normalized_deform_subset_intensities);

				// 归一化的"参考子区"和"变形子区"的灰度差
				const vec diff = normalized_refer_subset_intensities - normalized_deform_subset_intensities;

				// 变形子区点的灰度梯度
				deform_subset_gradient_x = deform_image_interpolator->GradientsX(deform_x, deform_y);
				deform_subset_gradient_y = deform_image_interpolator->GradientsY(deform_x, deform_y);

				// 伪逆矩阵
				PseudoInverseMatrix(deform_subset_gradient_x, deform_subset_gradient_y, delta_x, delta_y, dic_parameters.weight(), dic_parameters.shape_function_order(), pseudo);

				// 迭代增量
				dp = gg * pseudo * diff;
				warp_function = warp_function + dp;

				double u = dp(0);
				double v = dp(1);

				// 迭代次数加1
				++iter;
			}

			// 最终的相关系数
			ShapeFunction(x(r,c), y(r,c), delta_x, delta_y, warp_function, dic_parameters.shape_function_order(), deform_x, deform_y);

			//// 判断是否越界
			//if (deform_x.min() < 10 || deform_x.max() > 490 || deform_y.min() < 10 || deform_y.max() > 490)
			//{
			//	dic_output->set_warp_function(r, c, warp_function);
			//	dic_output->set_iteration_times(r, c, iter);
			//	dic_output->set_zncc(r, c, -1);
			//	dic_output->set_valid_sign(r, c, DICOutput::POI_FAIL);
			//	break;
			//}

			deform_subset_intensities = deform_image_interpolator->Values(deform_x, deform_y);

			const double zncc = ZNCC(refer_subset_intensities, deform_subset_intensities);

			// 记录计算结果
			dic_output->set_warp_function(r, c, warp_function);
			dic_output->set_iteration_times(r, c, iter);
			dic_output->set_zncc(r, c, zncc);

			if (zncc >= dic_parameters.zncc_threshold() && iter != dic_parameters.max_iteration_times())
				dic_output->set_valid_sign(r, c, DICOutput::POI_SUCCEED);
			else
				dic_output->set_valid_sign(r, c, DICOutput::POI_FAIL);

		}
	}

	delete deform_image_interpolator;

	return true;

}



bool RegisterSubpixelDisplacementWithBoundaryCheck(const arma::mat& refer_image, const arma::mat& deform_image, const DICParameters& dic_parameters, DICOutput* dic_output)
{
	// 网格点
	const mat x = dic_output->x();
	const mat y = dic_output->y();

	// ROI范围
	const QRect roi = dic_parameters.roi();

	mat mask(refer_image.n_rows, refer_image.n_cols, fill::zeros);
	for (int r = 0; r < mask.n_rows; ++r)
	{
		for (int c = 0; c < mask.n_cols; ++c)
		{
			if (r >= roi.top() && 
				r <= roi.bottom() && 
				c >= roi.left() && 
				c <= roi.right())
				mask(r, c) = 1;
		}
	}

	// "子区点"相对于"子区中心"的坐标
	mat delta_x, delta_y;
	SubsetPointsRelativeCoordinate(dic_parameters.subset_size(), delta_x, delta_y);

	// 变形图插值系数
	Interpolator* deform_image_interpolator = ConstructInterpolator(deform_image, dic_parameters.bspline_interpolation_order());

	int num_calculated = 0;
	const int tick = x.n_elem / 10;

#	pragma omp parallel for
	for (int r = 0; r < x.n_rows; ++r)
	{
		for (int c = 0; c < x.n_cols; ++c)
		{
			// 显示计算比例
#			pragma omp critical
			++num_calculated;
			if (num_calculated % tick == 0)
				printf("%.2f%%\t", 100. * num_calculated / x.n_elem);

			// 子区标识位
			const mat weight = SubsetIntensities(mask, x(r, c), y(r, c), dic_parameters.half_subset_size());

			// 子区标识位向量化
			const vec w = vectorise(weight);

			// 参考子区灰度
			const mat refer_subset_intensities = SubsetIntensities(refer_image, x(r, c), y(r, c), dic_parameters.half_subset_size());

			// 参考子区灰度向量化
			vec nr;
			NormalizeVectorize(refer_subset_intensities, weight, nr);

			// 迭代初值
			vec warp_function(dic_output->warp_function(r, c));

			// 迭代次数
			uword iter = 0;

			// 迭代
			vec dp = ones<vec>(warp_function.n_elem);	// 迭代增量

			while (!isConverge(dp, dic_parameters.error_threshold()) && iter <= dic_parameters.max_iteration_times())
			{
				// 确定变形子区点的"图像坐标"
				mat deform_x, deform_y;
				ShapeFunction(x(r, c), y(r, c), delta_x, delta_y, warp_function, dic_parameters.shape_function_order(), deform_x, deform_y);

				// 变形子区灰度
				const mat deform_subset_intensities = deform_image_interpolator->Values(deform_x, deform_y);
				
				// 变形子区灰度向量化
				vec nd;
				const double gg = NormalizeVectorize(deform_subset_intensities, weight, nd);

				// 归一化的"参考子区"和"变形子区"的灰度差
				const vec diff = nr - nd;

				// 变形子区点的灰度梯度
				const mat deform_subset_gradient_x = deform_image_interpolator->GradientsX(deform_x, deform_y);
				const mat deform_subset_gradient_y = deform_image_interpolator->GradientsY(deform_x, deform_y);

				// 伪逆矩阵
				mat pseudo;
				PseudoInverseMatrix(deform_subset_gradient_x, deform_subset_gradient_y, delta_x, delta_y, w, dic_parameters.shape_function_order(), pseudo);

				// 迭代增量
				dp = gg * pseudo * diff;
				warp_function = warp_function + dp;

				double u = dp(0);
				double v = dp(1);

				// 迭代次数加1
				++iter;
			}

			// 最终的相关系数
			mat deform_x, deform_y;
			ShapeFunction(x(r, c), y(r, c), delta_x, delta_y, warp_function, dic_parameters.shape_function_order(), deform_x, deform_y);
			const mat deform_subset_intensities = deform_image_interpolator->Values(deform_x, deform_y);
			const double zncc = ZNCC(refer_subset_intensities, deform_subset_intensities, weight);

			// 记录计算结果
			dic_output->set_warp_function(r, c, warp_function);
			dic_output->set_iteration_times(r, c, iter);
			dic_output->set_zncc(r, c, zncc);

			if (zncc >= dic_parameters.zncc_threshold() && iter != dic_parameters.max_iteration_times())
				dic_output->set_valid_sign(r, c, DICOutput::POI_SUCCEED);
			else
				dic_output->set_valid_sign(r, c, DICOutput::POI_FAIL);

		}
	}

	delete deform_image_interpolator;

	return true;



}



Interpolator* ConstructInterpolator(const mat& image, const int& n)
{
	Interpolator* interpolator;

	if (n == 3)
	{
		interpolator = new BicubicMOMSInterpolator(image);
	}
	else if (n == 5)
	{
		interpolator = new BiquinticBSplineInterpolatror(image);
	}
	else if (n == 7)
	{
		interpolator = new BisepticBSplineInterpolator(image);
	} 
	else
	{
		interpolator = new BiquinticBSplineInterpolatror(image);
	}

	return interpolator;
}



void SubsetPointsRelativeCoordinate(const int& m, mat& dx, mat& dy)
{
	dx.zeros(m, m);
	dy.zeros(m, m);

	const int n = (m - 1) / 2;
	for (int r = 0; r < m; ++r)
	{
		for (int c = 0; c < m; ++c)
		{
			dx(r, c) = c - n;
			dy(r, c) = r - n;
		}
	}

}



mat SubsetIntensities(const mat& image, const int& x0, const int& y0, const int& m) 
{ 
	return image(span(y0 - m, y0 + m), span(x0 - m, x0 + m)); 
}



bool isConverge(const vec& dp, const double& threshold)
{
	// 迭代增量
	const double dx = dp(0);
	const double dy = dp(1);
	const double dr = abs(dx) + abs(dy);

	// 收敛判据
	bool is_converge = false;
	if (dr < threshold)
		is_converge = true;

	return is_converge;
}



/*******************************************************************************************************************/
/**********************************            辅助计算函数              ********************************************/
/******************************************************************************************************************/



int ParameterTotalNumber(const int& order)
{
	int n = 0;

	if (order == 0)
		n = 2;
	else if (order == 1)
		n = 6;
	else if (order == 2)
		n = 12;
	//else if (order == 3)
	//	n = 20;
	//else if (order == 4)
	//	n = 30;
	//else if (order == 5)
	//	n = 42;

	return n;
}



double NormalizeVectorize(const mat& m, vec& v)
{
	// 向量化
	v = vectorise(m);

	// 零均值
	v = v - mean(v);

	// 归一化
	const double vv = norm(v, 2);
	v = v / vv;

	return vv;
}



double NormalizeVectorize(const mat& m, const mat& w, vec& v)
{
	// 将无效像素置零
	mat a = m % w;

	// 零均值
	a = a - accu(a) / accu(w);

	// 将无效像素置零
	a = a % w;

	// 归一化
	const double vv = sqrt(accu(a % a));
	a = a / vv;

	// 向量化
	v = vectorise(a);

	return vv;
}



double ZNCC(const mat& a, const mat& b)
{
	// 向量化
	vec aa = vectorise(a);
	vec bb = vectorise(b);

	// 零均值
	aa = aa - mean(aa);
	bb = bb - mean(bb);

	// 归一化
	aa = normalise(aa);
	bb = normalise(bb);

	// 相关系数
	const double zncc = dot(aa, bb);

	return zncc;
}



double ZNCC(const mat& a, const mat& b, const mat& w)
{
	vec aa;
	NormalizeVectorize(a, w, aa);

	vec bb;
	NormalizeVectorize(b, w, bb);

	const double zncc = dot(aa, bb);

	return zncc;
}