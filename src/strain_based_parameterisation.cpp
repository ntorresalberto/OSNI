#include "strain_based_parameterisation.h"




void StrainBasedParameterisation::updateParameterisation(const Eigen::VectorXd& t_qe,
                            const Eigen::VectorXd& t_dot_qe,
                            const Eigen::VectorXd& t_ddot_qe)
{
    Eigen::MatrixXd Phi(m_na, m_ne*m_na);

    for(unsigned int i=0; i<m_number_of_Chebyshev_points; i++){
        Phi = getPhi(m_chebyshev_points[i]);
        m_K_stack[i] = Phi*t_qe;
        m_dot_K_stack[i] = Phi*t_dot_qe;
        m_ddot_K_stack[i] = Phi*t_ddot_qe;
    }

    for(unsigned int i=0; i<m_number_of_Chebyshev_points; i++){
        m_Lambda_stack[i] = Eigen::Vector3d(1, 0, 0);
        m_dot_Lambda_stack[i] = Eigen::Vector3d::Zero();
        m_dot_Lambda_stack[i] = Eigen::Vector3d::Zero();
    }
}


Eigen::MatrixXd StrainBasedParameterisation::getPhi(const double& t_X,
                                                    const double& t_begin,
                                                    const double& t_end) const
{


    //  The coordinate must be transposed into the Chebyshev domain [-1, 1];
    double x = ( 2 * t_X - ( t_end + t_begin) ) / ( t_end - t_begin );

    //  Compute the values of the polynomial for every element of the strain field
    Eigen::VectorXd Phi_i(m_ne, 1);
    for(unsigned int i=0; i<m_ne; i++)
        Phi_i[i] = m_polynomial_base(i, x);


    //  Define the matrix of bases
    Eigen::MatrixXd Phi = Eigen::KroneckerProduct(Eigen::MatrixXd::Identity(m_na, m_na), Phi_i.transpose());


    return Phi;
}
