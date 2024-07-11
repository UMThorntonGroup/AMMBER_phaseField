// =================================================================================
// NUCLEATION FUNCTIONS
// =================================================================================

// =================================================================================
// Nucleation probability
// =================================================================================
template <int dim, int degree>
double customPDE<dim,degree>::getNucleationProbability(variableValueContainer variable_value, double dV, dealii::Point<dim> p, unsigned int variable_index) const
{
    unsigned int nuc_index = phase_index[variable_index];
    double G_nuc = fWell[nuc_index];
    double sum_op_sq = 0.0;
    for (unsigned int op_index = 0; op_index<num_ops; op_index++){
        sum_op_sq += variable_value(op_index)*variable_value(op_index);
    }
    std::vector<double> c(num_muFields, 0.0);
    std::vector<double> c_nuc(num_muFields, 0.0);
    // Calculate adjusted parabola minimum
    for (unsigned int mu_index=0; mu_index<num_muFields; mu_index++){
        double mu = variable_value(num_ops+mu_index);
        double mu_max = kWell[nuc_index][mu_index] * (1.0-cmin[nuc_index][mu_index]);
        double mu_min = kWell[nuc_index][mu_index] * (0.0-cmin[nuc_index][mu_index]);
        //std::cout << "mu:" << mu << " ";
        //std::cout << "mumin:" << mu_min << " ";
        //std::cout << "mumax:" << mu_max << " ";
        // mu = std::max(std::min(mu, mu_max), mu_min);
        double mu_nuc = std::max(std::min(mu, mu_max), mu_min); //note to future xander, some need mu, some need this.
        //std::cout << "munuc:" << mu_nuc << " ";
        c_nuc[mu_index] = cmin[nuc_index][mu_index]+mu_nuc/kWell[nuc_index][mu_index];
        //std::cout << "cnuc:" << c_nuc[mu_index] << " ";
        //c_nuc[mu_index] = std::max(std::min(c_nuc[mu_index], 1.0), 0.0);
        for (unsigned int op_index = 0; op_index<num_ops; op_index++){
            c[mu_index] += (variable_value(op_index)*variable_value(op_index)/sum_op_sq) // h
                          *(cmin[phase_index[op_index]][mu_index] + mu/(Va*kWell[phase_index[op_index]][mu_index]));
            //std::cout << "cmin:" << cmin[phase_index[op_index]][mu_index] << " ";
            //std::cout << "mu:" << mu << " ";
        }
        G_nuc += mu_nuc*(c_nuc[mu_index]-cmin[nuc_index][mu_index])/2.0;
        //std::cout << "G: " << G_nuc << " ";
    }
    // std::cout << "G" << variable_index << ": " << G_nuc << " ";
    G_nuc = -G_nuc;
    for (unsigned int op_index = 0; op_index<num_ops; op_index++){
        double h = (variable_value(op_index)*variable_value(op_index)/sum_op_sq);
        double g = fWell[phase_index[op_index]];
        for(unsigned int mu_index=0; mu_index<num_muFields; mu_index++){
            double mu = variable_value(num_ops+mu_index);
            //std::cout << "k:" << kWell[phase_index[op_index]][mu_index] << " ";
            //std::cout << "c:" << c[mu_index] << " ";
            //std::cout << "cmin:" << cmin[phase_index[op_index]][mu_index] << " ";
            //std::cout << "vv:" << variable_value(op_index) << " ";
            //std::cout << "ss:" << sum_op_sq << " ";
            g += 0.5*kWell[phase_index[op_index]][mu_index]*(c[mu_index]-cmin[phase_index[op_index]][mu_index])*(c[mu_index]-cmin[phase_index[op_index]][mu_index]);
            g += (mu*(c_nuc[mu_index]-c[mu_index]));
        }
        G_nuc += h*g;
    }
    //G_nuc = -G_nuc;
    // std::cout << "c_nuc" << variable_index << ": " << c_nuc[0] << " ";
    // std::cout << "GF" << variable_index << ": " << G_nuc << " ";
	// Calculate the nucleation rate
	double J=k1*exp(-k2/(std::max(G_nuc,1.0e-6)))/* *exp(-tau/(this->currentTime))*/;
    // std::cout << "J" << variable_index << ": " << J << " ";
	double retProb=1.0-exp(-J*userInputs.dtValue*((double)userInputs.steps_between_nucleation_attempts)*dV);
    // std::cout << "RP" << variable_index << ": " << retProb << " ";
    //std::cout<< G_nuc <<" ";
    return retProb;
}