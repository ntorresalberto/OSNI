/**
 * @file ode_a.h
 * @author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * @brief This files contains the abstraction for an ODE of the form \f[ \textbf{A} \dot{\textbf{x}} = \textbf{0} \f]
 * @version 0.1
 * @date 2022-07-07
 *
 * @copyright Copyright (c) 2022
 *
 */


#ifndef ODEA_H
#define ODEA_H

#include "spectral_integration_base.h"

#include "strain_based_parameterisation.h"

#include <unsupported/Eigen/KroneckerProduct>

namespace OSNI {



class ODEA : public ODESolverInterface
{
public:

    ODEA(const unsigned int t_state_dimension,
         const std::shared_ptr<const StrainBasedParameterisation> t_strain_parameterisation) : ODESolverInterface(t_state_dimension,
                                                                                                                       t_strain_parameterisation) {}

    ODEA(const unsigned int t_state_dimension,
         const std::shared_ptr<const StrainBasedParameterisation> t_strain_parameterisation,
         const unsigned int t_number_of_Chebyshev_points) : ODESolverInterface(t_state_dimension,
                                                                                    t_strain_parameterisation,
                                                                                    t_number_of_Chebyshev_points) {}


    virtual void integrate(const Eigen::VectorXd& t_initial_state) final;

    //virtual void getStateAtPoint(const unsigned int t_point) final;

    inline virtual Eigen::MatrixXd getStack()const final {return m_states_stack;}


    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point)=0;


private:

    const Eigen::MatrixXd m_D_NN { Eigen::KroneckerProduct(Eigen::MatrixXd::Identity(m_state_dimension, m_state_dimension), m_Dn_NN) };
    const Eigen::MatrixXd m_D_IN { Eigen::KroneckerProduct(Eigen::MatrixXd::Identity(m_state_dimension, m_state_dimension), m_Dn_IN) };

    Eigen::MatrixXd m_A_NN { m_D_NN };

    Eigen::MatrixXd m_A_at_Chebyshev_point { Eigen::MatrixXd::Zero(m_state_dimension, m_state_dimension) };

    void updateA();

    Eigen::VectorXd m_ivp { Eigen::VectorXd(m_state_dimension*(m_number_of_Chebyshev_points-1)) };

    Eigen::VectorXd m_states_stack { Eigen::VectorXd(m_state_dimension*(m_number_of_Chebyshev_points-1)) };


};


}   //  namespace OSNI



#endif // ODEA_H
