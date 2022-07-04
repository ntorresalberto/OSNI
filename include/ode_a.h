/**
 * @file ode_a.h
 * @author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * @brief This files contains the abstraction for an ODE of the form \f[ \textbf{A} \dot{\textbf{x}} = \textbf{0} \f]
 * @version 0.1
 * @date 03-07-2022
 *
 * @copyright Copyright (c) 2022
 *
 */


#ifndef ODE_A_H
#define ODE_A_H

#include "ode_solver_interface.h"




namespace OSNI {



class ODEA : public ODESolverInterface
{
public:

    ODEA(const unsigned int t_state_dimension) : m_state_dimension(t_state_dimension){}

    ODEA(const unsigned int t_state_dimension,
         const unsigned int t_number_of_Chebyshev_points) :  m_state_dimension(t_state_dimension),
                                                             m_number_of_Chebyshev_points(t_number_of_Chebyshev_points) {}




    virtual void integrate(const Eigen::VectorXd& t_initial_state) final;

    virtual Eigen::VectorXd getStateAtPoint(const unsigned int t_point)const final;

    inline virtual Eigen::MatrixXd getStack()const final {return m_states_stack; }

    virtual Eigen::MatrixXd computeMatrixAtChebyshevPoint(const unsigned int t_point)=0;

private:
    void updateA();


private:

    const unsigned int m_state_dimension;

    const unsigned int m_number_of_Chebyshev_points { 16 };

    const Eigen::MatrixXd m_D_NN {[&](){
            //  Get the Chebyshev differentiation matrix
            const Eigen::MatrixXd Dn = getDn(m_number_of_Chebyshev_points);

            //  Extract the block that define mutual infualces of the unknown states
            const auto Dn_NN = Dn.block(0, 0, m_number_of_Chebyshev_points-1, m_number_of_Chebyshev_points-1);

            //  Make a block diagonal matrix with the Dn_NN matrices
            return Eigen::KroneckerProduct(Eigen::MatrixXd::Identity(m_state_dimension, m_state_dimension), Dn_NN);
    }() };


    const Eigen::MatrixXd m_D_IN {[&](){
            //  Get the Chebyshev differentiation matrix
            const Eigen::MatrixXd Dn = getDn(m_number_of_Chebyshev_points);

            //  Extract the block that define influence of initial condition into the unknown states
            const auto Dn_IN = Dn.block(0, m_number_of_Chebyshev_points-1, m_number_of_Chebyshev_points-1, 1);

            //  Make a block diagonal matrix with the Dn_NN matrices
            return Eigen::KroneckerProduct(Eigen::MatrixXd::Identity(m_state_dimension, m_state_dimension), Dn_IN);
    }() };

    Eigen::MatrixXd m_A_NN { m_D_NN };

    Eigen::MatrixXd m_A_at_Chebyshev_point { Eigen::MatrixXd::Zero(m_state_dimension, m_state_dimension) };

    Eigen::VectorXd m_ivp { Eigen::VectorXd(m_state_dimension*(m_number_of_Chebyshev_points-1)) };

    Eigen::VectorXd m_states_stack { Eigen::VectorXd(m_state_dimension*(m_number_of_Chebyshev_points-1)) };


};


}   //  namespace OSNI



#endif // ODEA_H
