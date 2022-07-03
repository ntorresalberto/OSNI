#ifndef ODEAB_H
#define ODEAB_H

#include "spectral_integration_base.h"

#include "strain_based_parameterisation.h"

#include <unsupported/Eigen/KroneckerProduct>

namespace OSNI {


class ODEAb : public ODESolverInterface
{
public:
    ODEAb(const unsigned int t_state_dimension,
          const std::shared_ptr<const StrainBasedParameterisation> t_strain_parameterisation) : ODESolverInterface(t_state_dimension,
                                                                                                                       t_strain_parameterisation) {}

    ODEAb(const unsigned int t_state_dimension,
          const std::shared_ptr<const StrainBasedParameterisation> t_strain_parameterisation,
          const unsigned int t_number_of_Chebyshev_points) : ODESolverInterface(t_state_dimension,
                                                                                    t_strain_parameterisation,
                                                                                    t_number_of_Chebyshev_points) {}


    virtual void integrate(const Eigen::VectorXd& t_initial_state) final;

    //virtual void getStateAtPoint(const unsigned int t_point) final;

    inline virtual Eigen::MatrixXd getStack()const final {return m_states_stack;}


    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point)=0;

    virtual Eigen::VectorXd computerParametersVectorAtPoint(const unsigned int t_point)=0;


private:

    const Eigen::MatrixXd m_D_NN { Eigen::KroneckerProduct(Eigen::MatrixXd::Identity(m_state_dimension, m_state_dimension), m_Dn_NN) };
    const Eigen::MatrixXd m_D_IN { Eigen::KroneckerProduct(Eigen::MatrixXd::Identity(m_state_dimension, m_state_dimension), m_Dn_IN) };

    Eigen::MatrixXd m_A_NN { m_D_NN };
    Eigen::VectorXd m_b_NN { Eigen::VectorXd(m_state_dimension*(m_number_of_Chebyshev_points-1)) };


    Eigen::MatrixXd m_A_at_Chebyshev_point { Eigen::MatrixXd::Zero(m_state_dimension, m_state_dimension) };

    Eigen::VectorXd m_b_at_Chebyshev_point { Eigen::VectorXd::Zero(m_state_dimension) };


    void updateA();
    void updateb();

    Eigen::VectorXd m_ivp { Eigen::VectorXd(m_state_dimension*(m_number_of_Chebyshev_points-1)) };

    Eigen::VectorXd m_states_stack { Eigen::VectorXd(m_state_dimension*(m_number_of_Chebyshev_points-1)) };


};


}   //  namespace OSNI
#endif // ODEAB_H
