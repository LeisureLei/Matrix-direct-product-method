#include "optimize.h"

int main(int argc,char** argv)
{	
	vector<pose> gpspose,ousterpose;
	loadfile(gpspose,ousterpose);
	
	for(size_t i=stoi(argv[1]);i<stoi(argv[2])-1;i++)
	{	
		pose Ai,Bi;
		pose gpsi = findnearestneigbor(gpspose,ousterpose[i]);
		pose gpsiplus1 = findnearestneigbor(gpspose,ousterpose[i+1]);
		Bi.q = ousterpose[i].q.matrix().transpose()*ousterpose[i+1].q.matrix();
		Bi.t = ousterpose[i].q.matrix().transpose()*(ousterpose[i+1].t-ousterpose[i].t);
		Ai.q = gpsi.q.matrix().transpose()*gpsiplus1.q.matrix();
		Ai.t = gpsi.q.matrix().transpose()*(gpsiplus1.t-gpsi.t);
		A.push_back(Ai);
		B.push_back(Bi);
	}

	Eigen::Matrix3d I3 = Eigen::Matrix3d::Identity();	
	const int n = A.size()*9;
	Eigen::MatrixXd product(n,9);  
	for(size_t i=0;i<A.size();i++)
	{	//注意：矩阵按列展开为向量
		Eigen::Matrix<double,9,9> coeff = kroneckerproduct(A[i].q.matrix(),I3)-kroneckerproduct(I3,B[i].q.matrix().transpose());
		for(size_t j=0;j<coeff.rows();j++)
		{
			for(size_t k=0;k<coeff.cols();k++)
			{
				product(i*9+j,k) = coeff(j,k);  //对product逐个元素赋值
			}
		}
	}

	//计算旋转矩阵
	Eigen::JacobiSVD<Eigen::MatrixXd> svd ( product, Eigen::ComputeFullU|Eigen::ComputeFullV );
  const Eigen::MatrixXd V = svd.matrixV();
	Eigen::Matrix3d V9;
	V9<<V(0,8),V(1,8),V(2,8),V(3,8),V(4,8),V(5,8),V(6,8),V(7,8),V(8,8);
	Eigen::JacobiSVD<Eigen::MatrixXd> svd2 ( V9, Eigen::ComputeFullU|Eigen::ComputeFullV );
	const Eigen::MatrixXd u = svd2.matrixU();
	const Eigen::MatrixXd v = svd2.matrixV();
	Eigen::Matrix3d R = u*v.transpose();
	if(R.determinant() < 0){
		cout<<R.determinant()<<endl;
		R = -R;
	}
	cout<<"u*v^T:"<<endl<<R<<endl;
	cout<<"singularvalues:"<<endl<<svd2.singularValues()<<endl;
	cout<<"R determinant:"<<R.determinant()<<endl;

	//计算平移向量
	int n2 = A.size()*3;
	Eigen::MatrixXd leftA(n2,3);
	for(size_t i=0;i<A.size();i++)
	{
		Eigen::Matrix3d lefta(A[i].q.matrix()-I3);
		for(size_t j=0;j<lefta.rows();j++)
		{
			for(size_t k=0;k<lefta.cols();k++)
			{
				leftA(i*3+j,k) = lefta(j,k);
			}
		}
	}
	
	Eigen::VectorXd rightB(n2);
	for(size_t i=0;i<B.size();i++)
	{
		Eigen::Vector3d rightb(R*B[i].t-A[i].t);
		for(size_t j=0;j<rightb.size();j++)
		{
			rightB[i*3+j] = rightb[j];
		}	
	}
	
	Eigen::MatrixXd t = (leftA.transpose()*leftA).inverse()*leftA.transpose()*rightB;

	cout<<"t:"<<endl<<t<<endl;
	
	//ceres优化
	Eigen::Quaterniond q(R);
	double q_coeffs[4] = {q.w(),q.x(),q.y(),q.z()};
  double t_coeffs[3] = {t(0,0),t(1,0),t(2,0)};
	ceres::Problem problem;
	for(size_t i=0;i<A.size();i++)
	{
		ceres::CostFunction * costfunction = new ceres::AutoDiffCostFunction<CostFunction,6,4,3>(
		new CostFunction(A[i].q, B[i].q, A[i].t, B[i].t)
		);
		ceres::LossFunction * loss_function = new ceres::HuberLoss(1.0);
		problem.AddResidualBlock(costfunction, loss_function, q_coeffs, t_coeffs);
	}
	ceres::LocalParameterization* quaternionParameterization = new ceres::QuaternionParameterization;
  problem.SetParameterization(q_coeffs,quaternionParameterization);

  ceres::Solver::Options options;
  options.linear_solver_type = ceres::DENSE_SCHUR;
  options.max_num_iterations = 1000;
  ceres::Solver::Summary summary;
 	ceres::Solve(options, &problem, & summary);
  q = Eigen::Quaterniond(q_coeffs[0],q_coeffs[1],q_coeffs[2],q_coeffs[3]);
	cout<<q.matrix()<<endl;

	cout<<t_coeffs[0]<<" "<<t_coeffs[1]<<" "<<t_coeffs[2]<<endl;
	return 0;

}
