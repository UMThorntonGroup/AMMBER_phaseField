#include "../../include/matrixFreePDE.h"
#include <random>
#include <unordered_map>

template <int dim, int degree>
class customPDE: public MatrixFreePDE<dim,degree>
{
public:
    // Constructor
    customPDE(userInputParameters<dim> _userInputs): MatrixFreePDE<dim,degree>(_userInputs) , userInputs(_userInputs) {
        //Constructed variables & non-dimensionalization
        sigma *= 1.0/(l0*E0);
        l_gb *= 1.0/l0;
        Va *= 1.0/(l0*l0*l0);
        L *= E0*tau0;
        vecTimesX(D,tau0/(l0*l0));

        m0 = 6.0*sigma/l_gb;
        kappa = 0.75*sigma*l_gb;

        vecTimesX(fWell, 1.0/(l0*l0*l0*Va*E0));
        for(uint i=0;i<kWell.size(); ++i){
            vecTimesX(kWell[i], 1.0/(l0*l0*l0*Va*E0));
        }
        
        if(true){
            std::cout
                << "sigma: " << sigma << "\n" 
                << "l_gb: " << l_gb << "\n" 
                << "Va: " << Va << "\n" 
                << "L: " << L << "\n" 
                << "D[0]: " << D[0] << "\n" 
                << "m0: " << m0  << "\n" 
                << "kappa: " << kappa << "\n" 
                << "fWell[0]: " << fWell[0] <<  "\n"
                << "kWell[0][0]: " << kWell[0][ 0] << "\n";
        }
        
        //Defining seed
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        //Initializing distribution
        dist = distribution(0.0,1.0);
        //Initializing random variable
        rng = engine(seed);

    };

    // Function to set the initial conditions (in ICs_and_BCs.h)
    void setInitialCondition(const dealii::Point<dim> &p, const unsigned int index, double & scalar_IC, dealii::Vector<double> & vector_IC);

    // Function to set the non-uniform Dirichlet boundary conditions (in ICs_and_BCs.h)
    void setNonUniformDirichletBCs(const dealii::Point<dim> &p, const unsigned int index, const unsigned int direction, const double time, double & scalar_BC, dealii::Vector<double> & vector_BC);

    typedef std::mt19937_64 engine;
    typedef std::uniform_real_distribution<double> distribution;


private:
	#include "../../include/typeDefs.h"

	const userInputParameters<dim> userInputs;

	// Function to set the RHS of the governing equations for explicit time dependent equations (in equations.h)
    void explicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

    // Function to set the RHS of the governing equations for all other equations (in equations.h)
    void nonExplicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

	// Function to set the LHS of the governing equations (in equations.h)
	void equationLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

	// Function to set postprocessing expressions (in postprocess.h)
	#ifdef POSTPROCESS_FILE_EXISTS
	void postProcessedFields(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					variableContainer<dim,degree,dealii::VectorizedArray<double> > & pp_variable_list,
					const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;
	#endif
    // Function to set the nucleation probability (in nucleation.h)
    #ifdef NUCLEATION_FILE_EXISTS
    double getNucleationProbability(variableValueContainer variable_value, double dV) const;
    double getNucleationProbability(variableValueContainer variable_value, double dV, dealii::Point<dim> p, unsigned int variable_index) const;
    #endif
	
	// ================================================================
	// Methods specific to this subclass
	// ================================================================
    // Function to reshape a 1D vector to a 2D vector with n rows and m columns
    std::vector<std::vector<double>> reshapeVector(const std::vector<double>& inputVector, unsigned int n, unsigned int m) {
        int totalElements = n * m;

        // Check if the total number of elements matches the size of the input vector
        if (inputVector.size() != totalElements) {
            std::cerr << "[customPDE.h] Error: The total number of elements ["
                        << (n*m) << "] does not match the size of the input vector ["
                        << inputVector.size() << "]." << std::endl;
            std::exit(1);
            // Return an empty 2D vector to indicate an error
            return {};
        }

        // Reshape the vector
        std::vector<std::vector<double>> reshapedVector(n, std::vector<double>(m));

        for (int i = 0; i < totalElements; ++i) {
            reshapedVector[i / m][i % m] = inputVector[i];
        }

        return reshapedVector;
    }

    // key is phase ID, value is list of op indices
    std::unordered_map<unsigned int, std::vector<unsigned int> > make_map(std::vector<int> phase_id){
        std::unordered_map<unsigned int, std::vector<unsigned int>> out;
        for(unsigned int var_index = 0; var_index < phase_id.size(); ++var_index){
            out[phase_id[var_index]].push_back(var_index);
        }
        return out;
    }
    std::unordered_map<unsigned int, unsigned int> make_map_from_name(std::vector<std::string> phase_id, std::unordered_map<std::string, unsigned int> phase_index){
        std::unordered_map<unsigned int, unsigned int> out;
        for(unsigned int var_index = 0; var_index < phase_id.size(); ++var_index){
            out[phase_index[phase_id[var_index]]] = var_index;
        }
        return out;
    }
    std::vector<int> get_phase_index(unsigned int _num_phases, unsigned int _num_ops) const{
        std::vector<int> out(_num_ops);
        if(boost::iequals(userInputs.get_model_constant_string("set_ids_by"), "INDEX")){
            out = userInputs.get_model_constant_int_array("phase_id");
        }
        else if(boost::iequals(userInputs.get_model_constant_string("set_ids_by"), "NAME")){
            std::unordered_map<std::string, unsigned int> phase_id_by_name;
            std::vector<std::string> phase_names = userInputs.get_model_constant_string_array("phase_names");
            for(unsigned int i = 0; i<_num_phases; ++i){
                phase_id_by_name[phase_names[i]] = i;
            }
            std::vector<std::string> phase_id_names = userInputs.get_model_constant_string_array("phase_id");
            for(unsigned int op = 0; op<_num_ops; ++op){
                out[op] = phase_id_by_name[phase_id_names[op]];
            }
        }
        return out;
    }
    std::unordered_map<std::string, unsigned int> make_string_uint_map(std::vector<std::string> _phase_names, unsigned int _num_phases){
        std::unordered_map<std::string, unsigned int> phase_id_by_name;
        for(unsigned int i = 0; i<_num_phases; ++i){
            phase_id_by_name[_phase_names[i]] = i;
        }
        return phase_id_by_name;
    }
	// ================================================================
	// Model constants specific to this subclass
	// ================================================================
        unsigned int num_phases = userInputs.get_model_constant_int("num_phases");
        unsigned int num_comps  = userInputs.get_model_constant_int("num_comps");
        unsigned int num_muFields = num_comps-1;
        unsigned int num_ops = userInputs.get_model_constant_int("num_ops");
        std::vector<std::string> phase_names = userInputs.get_model_constant_string_array("phase_names");
        std::vector<int> phase_index = get_phase_index(num_phases, num_ops);

        // Values for non-dimensionalization
        double l0 = userInputs.get_model_constant_double("l0");
        double E0 = userInputs.get_model_constant_double("E0");
        double tau0 = userInputs.get_model_constant_double("tau0");
        
        // Phase Parameters
        std::vector<double> fWell = userInputs.get_model_constant_double_array("fWell");
        // 2D vectors of shape [num_phases][num_comps-1]
        std::vector<std::vector<double>> kWell = reshapeVector(userInputs.get_model_constant_double_array("kWell"),
                                                                  num_phases, num_muFields);
        std::vector<std::vector<double>> cmin  = reshapeVector(userInputs.get_model_constant_double_array("cmin"),
                                                                  num_phases, num_muFields);
 
        void vecTimesX(std::vector<double> &vec, double X){
          for(uint i=0; i<vec.size(); ++i){
            vec[i] *= X;
          }
        }
       
        // Model Parameters
        double sigma  = userInputs.get_model_constant_double("sigma");
        double l_gb   = userInputs.get_model_constant_double("l_gb");
        double Va     = userInputs.get_model_constant_double("Va");
        double L      = userInputs.get_model_constant_double("L");
        std::vector<double> D = userInputs.get_model_constant_double_array("D");

        double m0;
        double kappa;
        double gamma  = userInputs.get_model_constant_double("gamma");
        
        // Initial compositions if applicable
        std::vector<double> c0 = userInputs.get_model_constant_double_array("c0");
        
        // Declaring random number generator (Type std::mt19937_64)
        engine rng;
        // Declaring distribution (Type std::uniform_real_distribution<double>)
        distribution dist;
	// ================================================================

};
