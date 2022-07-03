#include "quaternion_integrator.h"

Eigen::MatrixXd QuaternionIntegrator::computeMatrixAtChebyshevPoint(const unsigned int t_point)
{
    m_K_at_point = m_strain_parameterisation->getKAtPoint(t_point);

    //  Compute the A matrix of Q' = 1/2 A(K) Q
    m_A_at_point  <<        0       ,   -m_K_at_point(0),   -m_K_at_point(1),   -m_K_at_point(2),
                    m_K_at_point(0) ,           0       ,    m_K_at_point(2),   -m_K_at_point(1),
                    m_K_at_point(1) ,   -m_K_at_point(2),           0       ,    m_K_at_point(0),
                    m_K_at_point(2) ,    m_K_at_point(1),   -m_K_at_point(0),           0       ;

    return 0.5*m_A_at_point;

}
