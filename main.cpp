#include <iostream>
#include <memory>

//#include "strain_based_parameterisation.h"

#include "ode_a.h"
#include "quaternion_integrator.h"

#include <boost/math/special_functions/chebyshev.hpp>


using namespace std;

int main()
{


//    const unsigned int number_of_Chebyshev_points = 16;

//    std::shared_ptr<StrainBasedParameterisation> strain_parameterisation = std::make_shared<StrainBasedParameterisation>(number_of_Chebyshev_points/*, [](const unsigned int t_point,
//                                                                                                                         const double t_x){return boost::math::chebyshev_t(t_point, t_x);}*/);

//    constexpr unsigned int ne = 3;
//    constexpr unsigned int na = 3;
//    Eigen::Matrix<double, ne*na, 1> qe;
//    Eigen::Matrix<double, ne*na, 1> dot_qe;
//    Eigen::Matrix<double, ne*na, 1> ddot_qe;
//    //  Here we give some value for the strain
//    qe <<   0,
//            0,
//            0,
//            1.2877691307032,
//           -1.63807499160786,
//            0.437406679142598,
//            0,
//            0,
//            0;

//    strain_parameterisation->updateParameterisation(qe, dot_qe, ddot_qe);

//    std::shared_ptr<OSNI::SpectralIntegrationBase> quaternion_integrator = std::make_shared<QuaternionIntegrator>(strain_parameterisation, number_of_Chebyshev_points);

//    Eigen::Vector4d q_init(1, 0, 0, 0);
//    quaternion_integrator->integrate(q_init);

//    std::cout << "Quaternion stack : \n" << quaternion_integrator->getStack();



    return 0;
}
