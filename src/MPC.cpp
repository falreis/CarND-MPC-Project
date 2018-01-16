#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"
#include "Helper.hpp"

using CppAD::AD;

//variables
size_t x_start = 0;                   //x = 0 + vt*cos(psi)*dt
size_t y_start = x_start + N;         //y = 0 + vt*sin(psi)*dt
size_t psi_start = y_start + N;       //psi = 0 + (v/Lf)*delta*dt
size_t v_start = psi_start + N;       //v = 0 + a*dt
size_t cte_start = v_start + N;       //cte = f(x)- y + (vt*sin(epsi)*dt)
size_t epsi_start = cte_start + N;    //epsi = psi - psi*des+ ((v/Lf)*delta*dt)
size_t delta_start = epsi_start + N;  //
size_t a_start = delta_start + N - 1; //
size_t ref_v = 10;                    //reference velocity (in m/s)

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.
    fg[0]= 0;

    const double weight_cte = 1.5;
    const double weight_epsi = 1;
    const double weight_velocity = 1;
    const double weight_delta = 5;
    const double weight_acc = 1;
    const double weight_gap_angle = 1;
    const double weight_gap_acc = 1;

    for (size_t t = 0; t < N; t++) {
      fg[0] += CppAD::pow(vars[cte_start + t], 2) * weight_cte;
      fg[0] += CppAD::pow(vars[epsi_start + t], 2) * weight_epsi;
      fg[0] += CppAD::pow(vars[v_start + t] - ref_v, 2) * weight_velocity;
    }

    // Minimize the use of actuators.
    for (size_t t = 0; t < N - 1; t++) {
      fg[0] += CppAD::pow(vars[delta_start + t], 2) * weight_delta;
      fg[0] += CppAD::pow(vars[a_start + t], 2) * weight_acc;
    }

    // Minimize the value gap between sequential actuations.
    for (size_t t = 0; t < N - 2; t++) {
      fg[0] += CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2) * weight_gap_angle;
      fg[0] += CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2) * weight_gap_acc;
    }

    fg[x_start+1] = vars[x_start];
    fg[y_start+1] = vars[y_start];
    fg[psi_start+1] = vars[psi_start];
    fg[v_start+1] = vars[v_start];
    fg[cte_start+1] = vars[cte_start];
    fg[epsi_start+1] = vars[epsi_start];

    for (size_t t = 1; t < N; t++) {
      // The state at time t+1 .
      AD<double> x1 = vars[x_start + t];
      AD<double> y1 = vars[y_start + t];
      AD<double> psi1 = vars[psi_start + t];
      AD<double> v1 = vars[v_start + t];
      AD<double> cte1 = vars[cte_start + t];
      AD<double> epsi1 = vars[epsi_start + t];

      // The state at time t.
      AD<double> x0 = vars[x_start + t - 1];
      AD<double> y0 = vars[y_start + t - 1];
      AD<double> psi0 = vars[psi_start + t - 1];
      AD<double> v0 = vars[v_start + t - 1];
      AD<double> cte0 = vars[cte_start + t - 1];
      AD<double> epsi0 = vars[epsi_start + t - 1];

      // Only consider the actuation at time t.
      AD<double> delta0 = vars[delta_start + t - 1];
      AD<double> a0 = vars[a_start + t - 1];

      //3rd order polynomial and its derivative derivative
      AD<double> f0 = (pow(x0,3) * coeffs[3]) + (pow(x0,2) * coeffs[2]) + (coeffs[1] * x0) + coeffs[0]; //ax^3 + bx^2 + cx + d
      AD<double> psides0 = CppAD::atan((3 * pow(x0,2) * coeffs[2]) + (2 * coeffs[1] * x0) + coeffs[0]); //3ax^2 + 2bx + c

      fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[1 + psi_start + t] = psi1 - (psi0 + v0 * delta0 / Lf * dt);
      fg[1 + v_start + t] = v1 - (v0 + a0 * dt);
      fg[1 + cte_start + t] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
      fg[1 + epsi_start + t] = epsi1 - ((psi0 - psides0) + v0 * delta0 / Lf * dt);
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() { }
MPC::~MPC() {}

void MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  //size_t i;   //unused
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  int n_vars = (N * 6) + ((N-1)*2);
  int n_constraints = (N * 6);

  Dvector vars(n_vars);
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);

  const double x = state[0];
  const double y = state[1];
  const double psi = state[2];
  const double v = state[3];
  const double cte = state[4];
  const double epsi = state[5];

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }
  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;

  // TODO: Set lower and upper limits for variables.
  for (size_t i = 0; i < delta_start; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  for (size_t i = delta_start; i < a_start; i++) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }

  for (int i = a_start; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  for (int j = 0; j < n_constraints; j++) {
    constraints_lowerbound[j] = 0;
    constraints_upperbound[j] = 0;
  }

  // object that computes objective and constraints
  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;

  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi; 

  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
        options
      , vars
      , vars_lowerbound
      , vars_upperbound
      , constraints_lowerbound
      , constraints_upperbound
      , fg_eval
      , solution
  );

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  this->steer_value = solution.x[delta_start];
  this->throttle_value = solution.x[a_start];

  this->x_vals.clear();
  this->y_vals.clear();

  const double xfactor = 8.0;

  for (int i = 0; i < N; ++i) {
    const double x_val = xfactor * i;
    const double y_val = (pow(x_val,3) * coeffs[3]) + (pow(x_val,2) * coeffs[2]) + (x_val * coeffs[1]) + coeffs[0];
    this->x_vals.push_back(x_val);
    this->y_vals.push_back(-y_val);
  }
}
