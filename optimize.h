#include <fstream>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/SVD>
#include <Eigen/Dense>
#include <ceres/ceres.h>
#include <sys/types.h>  //opendir
#include <dirent.h>  //readdir
#include <vector>
#include <regex>
#include <map>
#include <cmath>
#include <ctime>


using namespace std;



typedef struct{
	double time;
	Eigen::Quaterniond q;
	Eigen::Vector3d t;
}pose;

vector<pose> A,B;

void loadfile(vector<pose>& gpspose, vector<pose>& ousterpose)
{
	ifstream gpsfile("./gpspose.txt");
	if(!gpsfile.is_open())
	{
		cout<<"can't find gpspose.txt"<<endl;
		return ;
	}
	while(gpsfile)
	{	
		pose p;
		string s;
		getline(gpsfile,s);
		stringstream ss(s);
		ss>>p.time>>p.t[0]>>p.t[1]>>p.t[2]>>p.q.w()>>p.q.x()>>p.q.y()>>p.q.z();
		gpspose.push_back(p);
		//cout<<gpspose.size()<<":"<<p.time<<endl;
	}
	gpsfile.close();
	gpspose.pop_back();
	//cout<<"gpspose size:"<<gpspose.size()<<endl;

	ifstream ousterfile("./ousterpose.txt");
	if(!ousterfile.is_open())
	{
		cout<<"can't find ousterpose.txt"<<endl;
		return ;
	}
	while(ousterfile)
	{	
		pose q;
		string s;
		getline(ousterfile,s);
		stringstream ss(s);
		ss>>q.time>>q.t[0]>>q.t[1]>>q.t[2]>>q.q.w()>>q.q.x()>>q.q.y()>>q.q.z();
		ousterpose.push_back(q);
	}
	ousterfile.close();
	ousterpose.pop_back();
}

string gettimefromname(string name)
{
	regex rule("(\\d+)(.)(\\d+)");
	smatch sm;
	regex_search(name, sm, rule);
	return sm[0];
}

pose findnearestneigbor(const vector<pose>& gpspose,const pose& ousterposei)
{
	map<double,int> table;
	for(size_t i=0;i<gpspose.size();i++)
	{
		double dis = sqrt((gpspose[i].time - ousterposei.time)*(gpspose[i].time - ousterposei.time));
		table[dis] = i;
	}
	auto table_at = table.begin();
	pose nearestpose;
	nearestpose.time =  gpspose[table_at->second].time;
	nearestpose.t = gpspose[table_at->second].t;
	nearestpose.q = gpspose[table_at->second].q;
	return nearestpose;
}

Eigen::MatrixXd kroneckerproduct(const Eigen::MatrixXd& R1,const Eigen::MatrixXd& R2)
{
	Eigen::Matrix<double,9,9> product;
	product<<R1(0,0)*R2, R1(0,1)*R2 ,R1(0,2)*R2,
				R1(1,0)*R2, R1(1,1)*R2, R1(1,2)*R2,
					 R1(2,0)*R2, R1(2,1)*R2, R1(2,2)*R2;

	return product;
}

pose estimateK(const pose& TgpsiTogpsiplus1, const pose& TlidariTolidariplus1)
{
	Eigen::Matrix3d RA(TlidariTolidariplus1.q);
	Eigen::Matrix3d RB(TgpsiTogpsiplus1.q);
	Eigen::Vector3d tA(TlidariTolidariplus1.t);
	Eigen::Vector3d tB(TgpsiTogpsiplus1.t);
	Eigen::Matrix3d I3 = Eigen::Matrix3d::Identity();
	Eigen::MatrixXd A = kroneckerproduct(RA,I3)-kroneckerproduct(I3,RB.transpose());
	Eigen::EigenSolver<Eigen::MatrixXd> eigen_solver ( A.transpose()*A );
	cout << "matrix values = \n" << eigen_solver.eigenvalues() << endl;
  cout << "matrix vectors = \n" << eigen_solver.eigenvectors() << endl;
}

template<typename T>
void mat2RPY(const Eigen::Matrix<T, 3, 3>& m, T& roll, T& pitch, T& yaw)
{
    roll = atan2(m(2,1), m(2,2));
    pitch = atan2(-m(2,0), sqrt(m(2,1) * m(2,1) + m(2,2) * m(2,2)));
    yaw = atan2(m(1,0), m(0,0));
}

struct CostFunction
{
    private:
    Eigen::Vector3d _ta,_tb;
    Eigen::Quaterniond _Qa,_Qb;
    
    public:
    CostFunction(Eigen::Quaterniond Qa, Eigen::Quaterniond Qb,Eigen::Vector3d ta,Eigen::Vector3d tb ):
            _Qa(Qa), _Qb(Qb), _ta(ta), _tb(tb){}  

    template<typename T> 
    bool operator()(const T* const q4x1 , const T* const t3x1, T* residuals) const{
        Eigen::Quaternion<T> Rx( q4x1[0], q4x1[1],q4x1[2],q4x1[3]);
        Eigen::Matrix<T,3,1> tx;
        tx<<t3x1[0],t3x1[1],t3x1[2];

        Eigen::Quaternion<T> q_a = _Qa.cast<T>();  
        Eigen::Matrix<T,3,1>  t_a =  _ta.cast<T>();
        Eigen::Quaternion<T> q_b = _Qb.cast<T>();
        Eigen::Matrix<T, 3,1> t_b = _tb.cast<T>();

        Eigen::Matrix<T,3,3> R_odo = q_a.toRotationMatrix();
        Eigen::Matrix<T,3,3> Rrc = Rx.toRotationMatrix();

        Eigen::Matrix<T, 3,1> t_err = (R_odo - Eigen::Matrix<T,3,3>::Identity() ) * tx - (Rrc * t_b) + t_a;

        Eigen::Quaternion<T> q_err = Rx.conjugate() * q_a * Rx * q_b.conjugate();
        Eigen::Matrix<T,3,3> R_err = q_err.toRotationMatrix();

        T roll, pitch, yaw;
        mat2RPY(R_err,roll, pitch,yaw);

        residuals[0] = t_err[0];
        residuals[1] = t_err[1];
        residuals[2] = t_err[2];
        residuals[3] = T(5)*roll;
        residuals[4] = T(5)*pitch;
        residuals[5] = T(5)*yaw;

        return true;    
    }
};
