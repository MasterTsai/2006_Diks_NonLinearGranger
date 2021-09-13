#include <iostream>
#include <armadillo>
#include <cmath>
#include <Eigen/Dense>
using namespace std;
using namespace arma;
using namespace Eigen;

int main()
{
    vec a;
int h=5;

    MatrixXd b(5,5);
    MatrixXd im(5,5);

     a <<-2 << -1 << 0 << 1 << 2;

   for(int i1=0;i1<h;i1++)
   {

            b(i1,0)=1.0;
            b(i1,1)=a(i1);
            b(i1,2)=a(i1)*a(i1)+1.0/12.0;
            b(i1,3)=pow(a(i1),3)+a(i1)/4.0;
            b(i1,4)=1.0/80.0+pow(a(i1),2)/2.0+pow(a(i1),4);


    }

    im=b.inverse();

    cout<<im<<endl;
    return 0;
}