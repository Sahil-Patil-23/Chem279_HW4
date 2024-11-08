#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <sstream>
#include <fstream>
#include <armadillo>
#include <map>
using namespace std;


// Function that calculates the factorial of an input value- from HW2
int Factorial(int n){
    int factorial = 1;
    for(int i = 1; i <= n; i++){
        factorial *= i;
    }
    return factorial;
}


// Function that calculates the Double factorial of an input value- from HW2
int Double_Factorial(int n){
    int res = 1;

    for(int i = n; i > 1; i -= 2){
        res *= i;
    }

    return res;
}


// Function that calculates the binomial between 2 integers- from HW2
double Binomial(int n , int m){
    return Factorial(m) / ( Factorial(n) * Factorial(m - n) );
}


// **************************************************************************************************************************************************


// Calculates overlap in 1 dimension
double overlap_Integral_1D(double alpha, double beta, double center_a, double center_b, int lA, int lB) {
    // Calculate prefactor term
    double prefactor = exp(-alpha * beta * pow(center_a - center_b, 2) / (alpha + beta));
    prefactor *= sqrt(M_PI / (alpha + beta));

    // Calculate center
    double center_product = (alpha * center_a + beta * center_b) / (alpha + beta);

    // Summation of the angular momentum combinations
    double sum = 0.0;
    for (int i = 0; i <= lA; i++) {
        for (int j = 0; j <= lB; j++) {
            // Only (i + j) even terms contribute
            if ((i + j) % 2 == 0) {
                sum += Binomial(lA, i) * Binomial(lB, j) * (Double_Factorial(i + j - 1) 
                        * pow(center_product - center_a, lA - i) * pow(center_product - center_b, lB - j)) 
                        / pow(2 * (alpha + beta), double(i + j) / 2);
            }
        }
    }
    double integral = prefactor * sum;
    return integral;
}


// Calculates overlap in all 3 dimensions
double Overlap_Integral_3D(arma::rowvec centers_a, arma::rowvec centers_b, double alpha, double beta, arma::ivec lmn_a, arma::ivec lmn_b) {
    // Calculate the overlap integral in each dimension
    double integral = overlap_Integral_1D(alpha, beta, centers_a(0), centers_b(0), lmn_a(0), lmn_b(0)) *
                      overlap_Integral_1D(alpha, beta, centers_a(1), centers_b(1), lmn_a(1), lmn_b(1)) *
                      overlap_Integral_1D(alpha, beta, centers_a(2), centers_b(2), lmn_a(2), lmn_b(2));

    return integral;
}


// Struct to hold more intimate details on the molecules we are working with
struct Atomic_Orbitals{
    string AO_type;
    string A_symbol;
    int valence_elec;
    arma::rowvec center;
    arma::ivec L_M_N;
    arma::vec exponents;
    arma::vec contract_coefficients;
    arma::vec normalization_constants; 


    // Constructor Method
    Atomic_Orbitals(string AO_type, string A_symbol, int valence_elec, arma::rowvec center, arma::ivec L_M_N, arma::vec exponents, arma::vec contract_coefficients)
    : AO_type(AO_type), A_symbol(A_symbol), valence_elec(valence_elec), center(center), L_M_N(L_M_N), exponents(exponents), 
      contract_coefficients(contract_coefficients) {}

    // Calculate normalization constants- should be an upgrade to the way I did it in HW2
    void Calculate_Normalization_Constants(){
        double overlap1 = Overlap_Integral_3D(center, center, exponents(0), exponents(0), L_M_N, L_M_N);
        double overlap2 = Overlap_Integral_3D(center, center, exponents(1), exponents(1), L_M_N, L_M_N);
        double overlap3 = Overlap_Integral_3D(center, center, exponents(2), exponents(2), L_M_N, L_M_N);

        normalization_constants = {1.0/sqrt(overlap1) , 1.0/sqrt(overlap2) , 1.0/sqrt(overlap3)};
    }


};


// Function to calculate contracted overlap between 2 'Atomic_Orbitals' objects
double Calculate_Contracted_Overlap(Atomic_Orbitals A1, Atomic_Orbitals A2) {
    double contracted_overlap = 0.0;

    // Loop over all exponents
    for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
            double unnorm_overlap = Overlap_Integral_3D(A1.center, A2.center, 
                                                      A1.exponents(k), A2.exponents(l), 
                                                      A1.L_M_N, A2.L_M_N);
                                                      
            contracted_overlap += A1.contract_coefficients(k) * A2.contract_coefficients(l) * A1.normalization_constants(k) * A2.normalization_constants(l) * unnorm_overlap;
        }
    }
    return contracted_overlap;
}


// Structure to hold input file data
struct Atom{

    int _atomic_number;
    double _x, _y, _z;

};


// Function that will read in the file & create a vector of 'Atoms'
vector<Atom> Read_Atom_Para(string & filepath){

    ifstream file(filepath);
    if(!file){
        throw runtime_error("Unable to open file!"); // Program should let user know if file cant be accessed
    }

    vector<Atom> atoms; // Vector of class instances will be used throughout the program
    string line;

    getline(file, line);

    while(getline(file, line)){
        stringstream ss(line);
        
        Atom atom;
        
        ss >> atom._atomic_number >> atom._x >> atom._y >> atom._z; // Reading in all the needed info from the file
        atoms.push_back(atom);
    }

    return atoms;

}


// Function that calculates # of AO's
int Calculate_AO(vector<Atom> atoms){

    int AO = 0;

    for(auto atom : atoms){
        if(atom._atomic_number == 1){
            AO++;
        }
        else{
            AO += 4;
        }

    }

    return AO;
}


// Function that calculates # of valence electrons
int Calculate_valence_elec(vector<Atom> atoms){
    int valence = 0;

    for(auto atom : atoms){
        if(atom._atomic_number == 1){
            valence++;
        }
        else if(atom._atomic_number == 4){
            valence += 4;
        }
        else if(atom._atomic_number == 5){
            valence += 5;
        }
        else if(atom._atomic_number == 6){
            valence += 6;
        }
        else{
            valence += 7;
        }
    }

    return valence;
}


vector<Atomic_Orbitals> buildBasisFunctionsList(const vector<Atom>& atoms) {
    vector<Atomic_Orbitals> basisFunctionsList;

    // Define exponent and contraction maps
    map<string, arma::vec> exponentsMap = {
        {"H", {3.42525091, 0.62391373, 0.16885540}},
        {"C", {2.94124940, 0.68348310, 0.22228990}},
        {"N", {3.78045590, 0.87849660, 0.28571440}},
        {"O", {5.03315130, 1.16959610, 0.38038900}},
        {"F", {6.46480320, 1.50228120, 0.48858850}}
    };

    map<string, arma::vec> sContractionMap = {
        {"H", {0.15432897, 0.53532814, 0.44463454}},
        {"C", {-0.09996723, 0.39951283, 0.70011547}},
        {"N", {-0.09996723, 0.39951283, 0.70011547}},
        {"O", {-0.09996723, 0.39951283, 0.70011547}},
        {"F", {-0.09996723, 0.39951283, 0.70011547}}
    };

    map<string, arma::vec> pContractionMap = {
        {"C", {0.15591627, 0.60768372, 0.39195739}},
        {"N", {0.15591627, 0.60768372, 0.39195739}},
        {"O", {0.15591627, 0.60768372, 0.39195739}},
        {"F", {0.15591627, 0.60768372, 0.39195739}}
    };


 for (const auto& atom : atoms) {
        string atomSym;
        switch (atom._atomic_number) {
            case 1: atomSym = "H"; break;
            case 6: atomSym = "C"; break;
            case 7: atomSym = "N"; break;
            case 8: atomSym = "O"; break;
            case 9: atomSym = "F"; break;
            default: atomSym = "Unknown"; break;
        }

        int valence = Calculate_valence_elec({atom});
        arma::rowvec center = {atom._x, atom._y, atom._z};
        arma::ivec lmn_s = {0, 0, 0};

        // s-type function
        Atomic_Orbitals AO_s("1s", atomSym, valence, center, lmn_s, exponentsMap[atomSym], sContractionMap[atomSym]);
        AO_s.Calculate_Normalization_Constants();
        basisFunctionsList.push_back(AO_s);

        if (atomSym != "H") {
            // For elements like carbon, add p-type functions
            arma::ivec lmn_px = {1, 0, 0};
            arma::ivec lmn_py = {0, 1, 0};
            arma::ivec lmn_pz = {0, 0, 1};

            Atomic_Orbitals AO_2s("2s", atomSym, valence, center, lmn_s, exponentsMap[atomSym], pContractionMap[atomSym]);
            Atomic_Orbitals AO_px("px", atomSym, valence, center, lmn_px, exponentsMap[atomSym], pContractionMap[atomSym]);
            Atomic_Orbitals AO_py("py", atomSym, valence, center, lmn_py, exponentsMap[atomSym], pContractionMap[atomSym]);
            Atomic_Orbitals AO_pz("pz", atomSym, valence, center, lmn_pz, exponentsMap[atomSym], pContractionMap[atomSym]);

            AO_2s.Calculate_Normalization_Constants();
            AO_px.Calculate_Normalization_Constants();
            AO_py.Calculate_Normalization_Constants();
            AO_pz.Calculate_Normalization_Constants();

            basisFunctionsList.push_back(AO_2s);
            basisFunctionsList.push_back(AO_px);
            basisFunctionsList.push_back(AO_py);
            basisFunctionsList.push_back(AO_pz);
        }
    }

    return basisFunctionsList;
}


// *****************************************************************************************************************************************************


const double hartree_2eV_constant = 27.2114079527;

class CNDO{
    public:
        int num_AO;
        int num_valence_elec;
        vector<Atomic_Orbitals> basis_func;
        arma::mat Fock_Matrix;
        arma::mat Overlap_Matrix;
        arma::mat GammaMatrix;
        arma::mat hCoreMat;
        arma::mat AlphaFockMat;
        arma::mat BetaFockMat;
        arma::mat AlphaDensityMat;
        arma::mat BetaDensityMat;
        arma::mat OldAlphaDensityMat;
        arma::mat OldBetaDensityMat;
        arma::vec TotalDensity;
        arma::mat AlphaCoeffMat;
        arma::mat BetaCoeffMat;
        arma::vec AlphaEnergy;
        arma::vec BetaEnergy;
        arma::mat CoreHamil;
        map<int, int> aoIndexToAtom;
        map<string, map<string, double>> diagCNDOPara;
        map<string, double> offDiagCNDOPara;

        CNDO(vector<Atom>& atoms){

            basis_func = buildBasisFunctionsList(atoms);

            num_AO = Calculate_AO(atoms);
            num_valence_elec = Calculate_valence_elec(atoms);

            Fock_Matrix = arma::mat(num_AO, num_AO, arma::fill::zeros);
            Overlap_Matrix = arma::mat(num_AO, num_AO, arma::fill::zeros);
            AlphaCoeffMat = arma::mat(num_AO, num_AO, arma::fill::zeros);
            BetaCoeffMat = arma::mat(num_AO, num_AO, arma::fill::zeros);
            AlphaDensityMat = arma::mat(num_AO, num_AO, arma::fill::zeros);
            BetaDensityMat = arma::mat(num_AO, num_AO, arma::fill::zeros);
            TotalDensity = arma::vec(num_AO, arma::fill::zeros);
            AlphaEnergy = arma::vec(num_AO, arma::fill::zeros);
            BetaEnergy = arma::vec(num_AO, arma::fill::zeros);
            CoreHamil = arma::mat(num_AO, num_AO, arma::fill::zeros);

            // Initialize semi-empirical parameters
            const vector<string> atoms_CNDO = {"H", "C", "N", "O", "F"};
            const vector<string> Orbitals = {"1s", "2s", "px", "py", "pz"};
            const vector<double> Diag_Values = {7.176, 14.051, 5.572, 5.572, 5.572, 19.316, 7.275, 7.275, 7.275,
                                            25.390, 9.111, 9.111, 9.111, 32.272, 11.080, 11.080, 11.080};
            const vector<double> OffDiag_Values = {9., 21., 25., 31., 39.};


            // Map of semi-empirical parameters 
            // map<string, map<string, double>> diagCNDOPara;
            // map<string, double> offDiagCNDOPara;


            int index = 0;
            for (const string& atom : atoms_CNDO) {
                if (atom == "H") {
                    diagCNDOPara[atom]["1s"] = Diag_Values[index++];
                } else {
                    for (int i = 1; i < Orbitals.size(); ++i) {
                        diagCNDOPara[atom][Orbitals[i]] = Diag_Values[index++];
                    }
                }
                offDiagCNDOPara[atom] = OffDiag_Values[distance(atoms_CNDO.begin(), find(atoms_CNDO.begin(), atoms_CNDO.end(), atom))];
            }


            // Map the AO index to the atom it belongs to
            index = 0;
            for (int A = 0; A < atoms.size(); A++) {
                for (int i = 0; i < num_AO; i++) {
                    aoIndexToAtom[index++] = A;
                }
            }

        }


        // Calculate Overlap Matrix
        arma::mat Calculate_OverlapMatrix() {
            // Initialized the overlap matrix in constructor

            // Loop over all basis function combinations
            for (int i = 0; i < num_AO; i++) {
                for (int j = 0; j < num_AO; j++) {
                    Overlap_Matrix(i, j) = Calculate_Contracted_Overlap(basis_func[i], basis_func[j]);
                }
            }

            return Overlap_Matrix;
        }


        double pg2eIntegral(arma::rowvec center_a, arma::rowvec center_b, double sigmaA, double sigmaB) {
            double U = pow(M_PI * sigmaA, 1.5) * pow(M_PI * sigmaB, 1.5);
            double V2 = 1.0 / (sigmaA + sigmaB);

            double distance = norm(center_a - center_b, 2);

            if (distance == 0.0) {
                return U * sqrt(2 * V2) * sqrt(2 / M_PI) * hartree_2eV_constant;  // eq 3.15
            } 

            double sqrtT = sqrt(V2) * distance;
            double result = U / distance * erf(sqrtT); 
            return result * hartree_2eV_constant;
        }



        double Calculate_2E_Integral(Atomic_Orbitals A1, Atomic_Orbitals A2){
            if (!(accu(A1.L_M_N) == 0 && accu(A2.L_M_N) == 0)) {
                cout << "Error: 2e integrals is only for s orbitals" << endl;
                return 0.0;
            }
    
            // Compute the product of contraction coefficients and normalization constants
            arma::vec contractionNormProd1 = A1.contract_coefficients % A1.normalization_constants;
            arma::vec contractionNormProd2 = A2.contract_coefficients % A2.normalization_constants;

            int numExponents = A1.exponents.n_elem;
            double twoElectronIntegral = 0.0;

            // Compute the two-electron integral using nested loops
            for (int i = 0; i < numExponents; i++) {
                double exponent1_i = A1.exponents(i);
                double contractionNormProd1_i = contractionNormProd1(i);

                for (int j = 0; j < numExponents; j++) {
                    double exponent1_j = A1.exponents(j);
                    double contractionNormProd1_j = contractionNormProd1(j);
                    double sigma1 = 1.0 / (exponent1_i + exponent1_j);  // eq 3.10

                    for (int k = 0; k < numExponents; k++) {
                        double exponent2_k = A2.exponents(k);
                        double contractionNormProd2_k = contractionNormProd2(k);

                        for (int l = 0; l < numExponents; l++) {
                            double exponent2_l = A2.exponents(l);
                            double contractionNormProd2_l = contractionNormProd2(l);
                            double sigma2 = 1.0 / (exponent2_k + exponent2_l);  // eq 3.10

                            double primitiveIntegral = pg2eIntegral(A1.center, A2.center, sigma1, sigma2);  // eq 3.14
                            twoElectronIntegral += contractionNormProd1_i * contractionNormProd1_j * contractionNormProd2_k * contractionNormProd2_l * primitiveIntegral;
                        }
                    }
                }
            }

            return twoElectronIntegral;

        }


        // Function to calculate the Gamma Matrix
        arma::mat Calculate_GammaMatrix(){
            vector<Atomic_Orbitals> sBasisFunctionsList;

            // Loop through basis functions and select s orbitals
            for (const auto& basisFunc : basis_func) {
                // Check if the orbital is s-type
                if (basisFunc.L_M_N(0) == 0 && basisFunc.L_M_N(1) == 0 && basisFunc.L_M_N(2) == 0) { // Assuming lmn is a vector
                    sBasisFunctionsList.push_back(basisFunc);
                }
            }

            arma::mat gamma_matrix = arma::mat(num_AO, num_AO, arma::fill::zeros);

            // Loop over all s orbital basis function combinations
            for (size_t i = 0; i < sBasisFunctionsList.size(); i++) {
                for (size_t j = 0; j < sBasisFunctionsList.size(); j++) {
                    // Calculate the two-electron integral for the combination
                    gamma_matrix(i, j) = Calculate_2E_Integral(sBasisFunctionsList[i], sBasisFunctionsList[j]);
                }
            }

            return gamma_matrix;

        }


        // Calculate Total Density
        arma::vec Calculate_Total_Density(){

            for(int MU = 0; MU < num_AO; MU++){
                int A = aoIndexToAtom[MU];
                TotalDensity(A) += AlphaDensityMat(MU, MU) + BetaDensityMat(MU, MU);
            }

            return TotalDensity;
        }


        // Function that builds & calculates the Fock Matrix
        arma::mat Calculate_Fock_Matrix(arma::mat Total_Density, arma::mat gamma_mat, const vector<Atom>& atoms) {
        Fock_Matrix = arma::mat(num_AO, num_AO, arma::fill::zeros);
        // cout << "Fock_Matrix Size: " << Fock_Matrix.n_rows << " x " << Fock_Matrix.n_cols << endl;
        // cout << "Total_Density Size: " << Total_Density.n_rows << " x " << Total_Density.n_cols << endl;
        // cout << "gamma_mat Size: " << gamma_mat.n_rows << " x " << gamma_mat.n_cols << endl;
        
        // Loop over all AOs in the molecule
        for (int mu = 0; mu < num_AO; mu++) {
            for (int nu = 0; nu < num_AO; nu++) {
                int A = mu; 
                int B = nu; 

                string chemSymA = basis_func[A].A_symbol; // Accessing atomic symbol from basis function
                string chemSymB = basis_func[B].A_symbol; // Accessing atomic symbol from basis function

                double gammaAA = gamma_mat(A, A);
                double gammaAB = gamma_mat(A, B);
                double pAA = Total_Density(A);
                double ZA = Calculate_valence_elec({atoms[A]}); 

                // Calculate the diagonal elements of the matrix
                if (mu == nu) {
                    string AO_Type = basis_func[mu].AO_type; // Should be "1s" for hydrogen

                    // Check if the keys exist in the map before accessing them
                    if (diagCNDOPara.count(chemSymA) > 0 && diagCNDOPara[chemSymA].count(AO_Type) > 0) {
                        Fock_Matrix(mu, nu) = -diagCNDOPara.at(chemSymA).at(AO_Type) + 
                                            ((pAA - ZA) - (Total_Density(mu) - 0.5)) * gammaAA;

                        // Update the diagonal elements of the matrix when A != B
                        for (int B_idx = 0; B_idx < atoms.size(); B_idx++) {
                            if (A != B_idx) {
                                double pBB = Total_Density(B_idx);
                                double ZB = Calculate_valence_elec({atoms[B_idx]}); // Calculate valence of atom B
                                Fock_Matrix(mu, nu) += (pBB - ZB) * gamma_mat(A, B_idx); // Update Fock matrix
                            }
                        }
                    } else {
                        // cerr << "Key not found in diagCNDOPara for " << chemSymA << " and AO type " << AO_Type << endl;
                    }
                } else { // Calculate the off-diagonal elements of the matrix
                    if (offDiagCNDOPara.count(chemSymA) > 0 && offDiagCNDOPara.count(chemSymB) > 0) {
                        Fock_Matrix(mu, nu) = (-offDiagCNDOPara.at(chemSymA) - offDiagCNDOPara.at(chemSymB)) / 2.0 * Overlap_Matrix(mu, nu) - (Total_Density(mu) * gammaAB);
                    } else {
                        // cerr << "Key not found in offDiagCNDOPara for " << chemSymA << " or " << chemSymB << endl;
                    }
                }
            }
        }
        return Fock_Matrix;
    }


    // Function that builds the Core Hamiltonian Matrix
    arma::mat Caclulate_CoreHamil_Matrix(const vector<Atom>& atoms, arma::mat gamma_mat){

        for(int mu = 0; mu < num_AO; mu++){
            for(int nu = 0; nu < num_AO; nu++){
                int A = mu; 
                int B = nu; 

                // Accessing atomic symbol from basis function
                string chemSymA = basis_func[A].A_symbol;
                string chemSymB = basis_func[B].A_symbol;

                // 
                double gammaAA = gamma_mat(A, A);
                double gammaAB = gamma_mat(A, B);
                double ZA = Calculate_valence_elec({atoms[A]});

                if(mu == nu){
                    string AO_type = basis_func[mu].AO_type;
                    CoreHamil(mu, nu) = -diagCNDOPara.at(chemSymA).at(AO_type) - 
                                        (ZA - 0.5) * gammaAA;

                    for(int B = 0; B < atoms.size(); B++){
                        if(A != B){
                            double ZB = Calculate_valence_elec({atoms[B]});
                            double gammAB = gamma_mat(A, B);
                            CoreHamil(mu, nu) -= ZB * gammaAB; 
                        }
                    }
                }
                else{
                    CoreHamil(mu, nu) = (-offDiagCNDOPara.at(chemSymA) -offDiagCNDOPara.at(chemSymB)) / (2.0 * Overlap_Matrix(mu, nu));
                }
            }
        }
        return CoreHamil;
    }


    // Function that builds and calculates the Density Matrix
    arma::mat Calculate_Density_Matrix(arma::mat coeff_matrix_A, int type){
        // if type is 0 -> alpha
        // if type is 1 -> beta

        arma::mat Density_Matrix;
        Density_Matrix = arma::mat(num_AO, num_AO, arma::fill::zeros);

        if(type == 0){
            for(int mu = 0; mu < num_AO; mu++){
                for(int nu = 0; nu < num_AO; nu++){
                    for(int i = 0; i < (num_valence_elec / (2 + (num_valence_elec % 2))); i++){
                        Density_Matrix(mu, nu) += AlphaCoeffMat(mu, i) * BetaCoeffMat(nu, i);
                    }
                }
            }
        }
        else if(type == 1){
            for(int mu = 0; mu < num_AO; mu++){
                for(int nu = 0; nu < num_AO; nu++){
                    for(int i = 0; i < (num_valence_elec / 2); i++){
                        Density_Matrix(mu, nu) += AlphaCoeffMat(mu, i) * BetaCoeffMat(nu, i);
                    }
                }
            }
        }
        return Density_Matrix;
    }


    // Function that calculates the total energy
    double Calculate_Total_Energy(arma::mat AlphaDensityMat, arma::mat CoreHamil, arma::mat AlphaFockMat, arma::mat BetaFockMat, arma::mat BetaDensityMat){
        double total_E = 0.0;
        for(int mu = 0; mu < num_AO; mu++){
            for(int nu = 0; nu < num_AO; nu++){
                total_E = (AlphaDensityMat(mu, nu) * (CoreHamil(mu,nu) + AlphaFockMat(mu,nu)) + BetaDensityMat(mu,nu) * (CoreHamil(mu,nu) + BetaFockMat(mu,nu)));
            }
        }
        total_E /= 2.0;
        return total_E;
    }


    // Function that runs the SCF Cycle & updates relevant matrices
    void SCF_Cycle(arma::mat gamma_mat, const vector<Atom>& atoms){

        bool converge = false;
        int SCF_cycle_count = 0;

        while(!converge){
            SCF_cycle_count++;

            AlphaFockMat = Calculate_Fock_Matrix(AlphaDensityMat, gamma_mat, atoms);
            BetaFockMat = Calculate_Fock_Matrix(BetaDensityMat, gamma_mat, atoms);

            eig_sym(AlphaEnergy, AlphaCoeffMat, AlphaFockMat);
            eig_sym(BetaEnergy, BetaCoeffMat, BetaFockMat);

            OldAlphaDensityMat = AlphaDensityMat;
            OldBetaDensityMat = BetaDensityMat;

            AlphaDensityMat = Calculate_Density_Matrix(AlphaCoeffMat, 0);
            BetaDensityMat = Calculate_Density_Matrix(BetaCoeffMat, 1);
            TotalDensity = Calculate_Total_Density();


            // Print matrices
            cout << "Alpha Density Matrix: " << endl;
            cout << AlphaDensityMat << endl;
            cout << "Beta Density Matrix: " << endl;
            cout << BetaDensityMat << endl;
            cout << "Total Density: " << endl;
            cout << TotalDensity << endl;
            cout << "Alpha Fock Matrix: " << endl;
            cout << AlphaFockMat << endl;
            cout << "Beta Fock Matrix: " << endl;
            cout << BetaFockMat << endl;
            cout << "Alpha Coeff Matrix: " << endl;
            cout << AlphaCoeffMat << endl;
            cout << "Beta Coeff Matrix: " << endl;
            cout << BetaCoeffMat << endl;

            if(abs(AlphaDensityMat - OldAlphaDensityMat).max() < 1e-6 && abs(BetaDensityMat - OldBetaDensityMat).max() < 1e-6){
                converge = true;

                cout << "*************************************************" << endl;
                cout << "SCF cycle converged after " << SCF_cycle_count << " iterations!" << endl;
                double energy = Calculate_Total_Energy(AlphaDensityMat, CoreHamil, AlphaFockMat, BetaFockMat, BetaDensityMat);
                cout << "Total Energy: " << energy << " eV" << endl;



            }
        }

    }


};



int main(){

    // Using 'H2.txt'
    string filepath_H2 = "../sample_input/H2.txt";
    vector<Atom> H2 = Read_Atom_Para(filepath_H2);
    CNDO H2_CNDO = CNDO(H2);
    arma::mat H2_overlapmat = H2_CNDO.Calculate_OverlapMatrix();
    arma::mat H2_Gammamat = H2_CNDO.Calculate_GammaMatrix();
    arma::mat H2_Fock = H2_CNDO.Calculate_Fock_Matrix(H2_CNDO.TotalDensity, H2_Gammamat, H2);
    arma::mat H2_CoreHamil = H2_CNDO.Caclulate_CoreHamil_Matrix(H2, H2_Gammamat);
    cout << "gamma" << endl;
    cout << H2_Gammamat << endl;
    cout << "Overlap" << endl;
    cout << H2_overlapmat << endl;
    cout << "H_core" << endl;
    cout << H2_Fock << endl;

    H2_CNDO.SCF_Cycle(H2_Gammamat, H2);

    cout << endl << endl;


    // Using 'HF.txt'
    string filepath_HF = "../sample_input/HF.txt";
    vector<Atom> HF = Read_Atom_Para(filepath_HF);
    CNDO HF_CNDO = CNDO(HF);
    arma::mat HF_overlapmat = HF_CNDO.Calculate_OverlapMatrix();
    arma::mat HF_Gammamat = HF_CNDO.Calculate_GammaMatrix();
    arma::mat HF_Fock = HF_CNDO.Calculate_Fock_Matrix(HF_CNDO.TotalDensity, HF_Gammamat, HF);
    cout << "gamma" << endl;
    cout << HF_Gammamat[0] << "  " << HF_Gammamat[1] << endl;
    cout << HF_Gammamat[1] << "  " << HF_Gammamat[0] << endl;
    cout << endl;

    cout << "Overlap" << endl;
    cout << HF_overlapmat << endl;
    cout << "H_Core" << endl;
    cout << HF_Fock << endl;

    HF_CNDO.SCF_Cycle(HF_Gammamat, HF);

    cout << endl << endl;

   



    // Using 'HO.txt'
    string filepath_HO = "../sample_input/HO.txt";
    vector<Atom> HO = Read_Atom_Para(filepath_HO);
    CNDO HO_CNDO = CNDO(HO);
    arma::mat HO_overlapmat = HO_CNDO.Calculate_OverlapMatrix();
    arma::mat HO_Gammamat = HO_CNDO.Calculate_GammaMatrix();
    arma::mat HO_Fock = HO_CNDO.Calculate_Fock_Matrix(HO_CNDO.TotalDensity, HO_Gammamat, HO);
    cout << "gamma" << endl;
    cout << HO_Gammamat[0] << "  " << HO_Gammamat[1] << endl;
    cout << HO_Gammamat[1] << "  " << HO_Gammamat[0] << endl;
    cout << endl;
    cout << "Overlap" << endl;
    cout << HO_overlapmat << endl;
    cout << "H_Core" << endl;
    cout << HO_Fock << endl;

    HO_CNDO.SCF_Cycle(HO_Gammamat, HO);

    cout << endl << endl;
    
    return 0;
}