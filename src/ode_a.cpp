#include "ode_a.h"

//#include "chebyshev_differentiation.hpp"

//#include <unsupported/Eigen/KroneckerProduct>

namespace OSNI {



void ODEA::integrate(const Eigen::VectorXd& t_initial_state)
{
    updateA();

    m_ivp = m_D_IN*t_initial_state;

    m_states_stack = - m_A_NN.inverse() * m_ivp;

}

Eigen::VectorXd ODEA::getStateAtPoint(const unsigned int t_point)const
{

    Eigen::VectorXd state_at_point(m_state_dimension);

    for(unsigned int i=0; i<m_state_dimension; i++){
        state_at_point(i) = m_states_stack(t_point + i*(m_number_of_Chebyshev_points-1));
    }

    return state_at_point;
}




void ODEA::updateA()
{
    for(unsigned int current_chebyshev_point=0; current_chebyshev_point<m_number_of_Chebyshev_points-1; current_chebyshev_point++){

        //  Compute the matrix at the Chebyshev point
        m_A_at_Chebyshev_point = computeMatrixAtChebyshevPoint(current_chebyshev_point);

        //  Perform the update of the A matrix
        int row_index, col_index;
        for (unsigned int row = 0; row < m_state_dimension; ++row) {
            for (unsigned int col = 0; col < m_state_dimension; ++col) {

                row_index = row*(m_number_of_Chebyshev_points-1) + current_chebyshev_point;
                col_index = col*(m_number_of_Chebyshev_points-1) + current_chebyshev_point;

                m_A_NN(row_index, col_index) = m_D_NN(row_index, col_index) - m_A_at_Chebyshev_point(row, col);
            }
        }

    }

}


}
