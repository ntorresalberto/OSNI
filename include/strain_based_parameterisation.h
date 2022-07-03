#ifndef STRAINBASEDPARAMETERISATION_H
#define STRAINBASEDPARAMETERISATION_H

#include <Eigen/Dense>
#include <vector>



#include <boost/math/special_functions/legendre.hpp>
#include <unsupported/Eigen/KroneckerProduct>


/*!
 * \brief ComputeChebyshevPoints Computes the Chebyshev points in the given interval
 * \tparam t_N The number of Chebyshev points.
 * \tparam t_L The length of the interval. Default 1 for the interval [0, 1]
 * \return An std::array containing the Chebyshev points
 */
static std::vector<double> ComputeChebyshevPoints(const unsigned int t_number_of_Chebyshev_points,
                                                                                const unsigned int t_L=1)
{
    std::vector<double> x(t_number_of_Chebyshev_points);

    unsigned int j = 0;
    std::generate(x.begin(), x.end(), [&](){
        return (static_cast<double>(t_L)/2)*(1 +cos( M_PI * static_cast<double>(j++) / static_cast<double>(t_number_of_Chebyshev_points-1) ));
    });

    return x;
}




class StrainBasedParameterisation
{
public:
    StrainBasedParameterisation()=default;

    StrainBasedParameterisation(const unsigned int t_number_of_Chebyshev_points) :
                m_number_of_Chebyshev_points(t_number_of_Chebyshev_points) {}

    StrainBasedParameterisation(const unsigned int t_number_of_Chebyshev_points,
                                std::function<double(const unsigned int, const double&)> t_polynomial_base) :
                            m_number_of_Chebyshev_points(t_number_of_Chebyshev_points),
                            m_polynomial_base(t_polynomial_base){}


    void updateParameterisation(const Eigen::VectorXd& t_qe,
                                const Eigen::VectorXd& t_dot_qe,
                                const Eigen::VectorXd& t_ddot_qe);

    inline Eigen::Vector3d getKAtPoint(const unsigned int t_point) const {return m_K_stack[t_point]; }

    inline Eigen::Vector3d getLambdaAtPoint(const unsigned int t_point) const {return m_Lambda_stack[t_point]; }

private:

    const unsigned int m_number_of_Chebyshev_points { 16 };

    const std::vector<double> m_chebyshev_points { ComputeChebyshevPoints(m_number_of_Chebyshev_points) };


    const unsigned int m_ne { 3 };
    const unsigned int m_na { 3 };

    const unsigned int m_n { m_ne*m_na };

    Eigen::VectorXd m_qe { Eigen::VectorXd::Zero(m_n) };
    Eigen::VectorXd m_dot_qe { Eigen::VectorXd::Zero(m_n) };
    Eigen::VectorXd m_ddot_qe { Eigen::VectorXd::Zero(m_n) };


    Eigen::MatrixXd getPhi(const double& t_X,
                            const double& t_begin=0,
                            const double& t_end=1) const;

    std::function<double(const unsigned int, const double&)> m_polynomial_base { [](const unsigned int t_point, const double& t_x) {return boost::math::legendre_p(t_point, t_x);} };


    std::vector<Eigen::Vector3d> m_K_stack { std::vector<Eigen::Vector3d>( m_number_of_Chebyshev_points ) };
    std::vector<Eigen::Vector3d> m_dot_K_stack { std::vector<Eigen::Vector3d>( m_number_of_Chebyshev_points ) };
    std::vector<Eigen::Vector3d> m_ddot_K_stack { std::vector<Eigen::Vector3d>( m_number_of_Chebyshev_points ) };

    std::vector<Eigen::Vector3d> m_Lambda_stack { std::vector<Eigen::Vector3d>( m_number_of_Chebyshev_points ) };
    std::vector<Eigen::Vector3d> m_dot_Lambda_stack { std::vector<Eigen::Vector3d>( m_number_of_Chebyshev_points ) };
    std::vector<Eigen::Vector3d> m_ddot_Lambda_stack { std::vector<Eigen::Vector3d>( m_number_of_Chebyshev_points ) };


};

#endif // STRAINBASEDPARAMETERISATION_H
