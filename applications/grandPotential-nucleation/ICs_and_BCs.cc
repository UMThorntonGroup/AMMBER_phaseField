// ===========================================================================
// FUNCTION FOR INITIAL CONDITIONS
// ===========================================================================

template <int dim, int degree>
void customPDE<dim,degree>::setInitialCondition(const dealii::Point<dim> &p, const unsigned int index, double & scalar_IC, dealii::Vector<double> & vector_IC){
    // ---------------------------------------------------------------------
    // ENTER THE INITIAL CONDITIONS HERE
    // ---------------------------------------------------------------------
    // Enter the function describing conditions for the fields at point "p".
    // Use "if" statements to set the initial condition for each variable
    // according to its variable index

    // Precalculating everything makes writing initial conditions easier. May take slightly more runtime.
    std::vector<double> op_vals(num_ops, 0.0);
    std::vector<double> mu_vals(num_muFields, 0.0);


    // double center[3] = {5.0,5.0,5.0};
    // double rad = 1.0;
    // double xwidth = 2.5;
    // double dist2 = 0.0;
    // double dist = 0.0;
    // double xdist = p[0]-center[0];
    // for(unsigned int xyz = 0; xyz<dim; xyz++){
    //     dist2 += (p[xyz]-center[xyz])*(p[xyz]-center[xyz]);
    // }
    // dist = std::sqrt(dist2);
    // double tanh_profile = 0.5*(1+std::tanh(-2.0*(dist-rad)));
    // double tanh_profileB= 0.5*(1+std::tanh(-2.0*(std::abs(xdist)-xwidth)));


    op_vals[0] = 1.0;//liquid
    //mu_vals[0] = 12.0;//component 0
    for(unsigned int mu_index = 0; mu_index<num_muFields; ++mu_index){
        mu_vals[0] = kWell[0][mu_index]*(c0[mu_index]-cmin[0][mu_index]);
    }

    // Submit fields
    scalar_IC = 0.0;
    for(unsigned int op_index = 0; op_index<num_ops; ++op_index){
        if(index==op_index){scalar_IC = op_vals[op_index];}
    }
    for(unsigned int mu_index = 0; mu_index<num_muFields; ++mu_index){
        if(index==num_ops+mu_index){scalar_IC = mu_vals[mu_index];}
    }
}

// ===========================================================================
// FUNCTION FOR NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS
// ===========================================================================

template <int dim, int degree>
void customPDE<dim,degree>::setNonUniformDirichletBCs(const dealii::Point<dim> &p, const unsigned int index, const unsigned int direction, const double time, double & scalar_BC, dealii::Vector<double> & vector_BC)
{
    // --------------------------------------------------------------------------
    // ENTER THE NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS HERE
    // --------------------------------------------------------------------------
    // Enter the function describing conditions for the fields at point "p".
    // Use "if" statements to set the boundary condition for each variable
    // according to its variable index. This function can be left blank if there
    // are no non-uniform Dirichlet boundary conditions. For BCs that change in
    // time, you can access the current time through the variable "time". The
    // boundary index can be accessed via the variable "direction", which starts
    // at zero and uses the same order as the BC specification in parameters.in
    // (i.e. left = 0, right = 1, bottom = 2, top = 3, front = 4, back = 5).


    // -------------------------------------------------------------------------

}
