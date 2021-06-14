#include <cmath>
#include <math.h>
#include <vector>
#include <numeric>

class Viscosity {
    private:
        double viscosityMixture = 0.0; 
        double fluidity = 0.0;
        double E = 0.0; 
        const double A = 0.375;     
        std::vector<double> phi;  
        std::vector<double> wilkeNumerator;
        double wilkeDenumerator = 0.0;
        std::vector<double> wilkeDenumeratorSum;
        double mixtureMolarMass = 0.0; 
        double mixtureFractionDenumerator = 0.0; 
        std::vector<double> moleFractions; 
        std::vector<double> mixtureFractions; 
    public: 
        std::vector<double> massToMoleFractions(std::vector<double> &massFractions, std::vector<double> &molarMasses, int numberOfComponents);
        double wilkeViscosity(std::vector<double> &massFractions, std::vector<double> &molarMasses, std::vector<double> &viscosities, int numberOfComponents);
        double davidsonViscosity(std::vector<double> &massFractions, std::vector<double> &molarMasses, std::vector<double> &viscosities, int numberOfComponents);
}; 

std::vector<double> Viscosity::massToMoleFractions(std::vector<double> &massFractions, std::vector<double> &molarMasses, int numberOfComponents){
    std::cout << "Mixture mole fractions: ";

    for(int i = 0; i < numberOfComponents; i++){
        mixtureMolarMass += massFractions[i] / molarMasses[i]; 
    }

    for(int i = 0; i < numberOfComponents; i++){
        moleFractions.push_back((massFractions[i] / molarMasses[i]) / mixtureMolarMass);
        std::cout << moleFractions[i] << " ";
    }

    std::cout << "\n";
    return moleFractions; 
}

double Viscosity::wilkeViscosity(std::vector<double> &massFractions, std::vector<double> &molarMasses, std::vector<double> &viscosities, int numberOfComponents){
    massToMoleFractions(massFractions, molarMasses, numberOfComponents); 
    
    for(int i = 0; i < numberOfComponents; i++){
        for(int j = 0; j < numberOfComponents; j++){
            phi.push_back(pow((1 + pow(viscosities[i] / viscosities[j],0.5) * pow(molarMasses[j] / molarMasses[i],0.25)),2) / pow(8 * (1 + molarMasses[i] / molarMasses[j]),0.5)); 
            wilkeDenumerator += moleFractions[j] * phi[j];   
        }
        wilkeDenumeratorSum.push_back(wilkeDenumerator); 
        wilkeDenumerator = 0.0; 
        phi.clear();
        wilkeNumerator.push_back(moleFractions[i] * viscosities[i]); 
        viscosityMixture += wilkeNumerator[i] / wilkeDenumeratorSum[i];  
    }
    return viscosityMixture; 
}

double Viscosity::davidsonViscosity(std::vector<double> &massFractions, std::vector<double> &molarMasses, std::vector<double> &viscosities, int numberOfComponents){
    massToMoleFractions(massFractions, molarMasses, numberOfComponents);

    for(int i = 0; i < numberOfComponents; i++){
        mixtureFractionDenumerator += moleFractions[i] * sqrt(molarMasses[i]);
    }

    for(int j = 0; j < numberOfComponents; j++){
        mixtureFractions.push_back((moleFractions[j] * sqrt(molarMasses[j])) / mixtureFractionDenumerator); 
    }
    
    for(int i = 0; i < numberOfComponents; i++){
        for(int j = 0; j < numberOfComponents; j++){
                E = (2*sqrt(molarMasses[i]) * sqrt(molarMasses[j])) / (molarMasses[i] + molarMasses[j]);
                fluidity += ((mixtureFractions[i] * mixtureFractions[j]) / (sqrt(viscosities[i]) * sqrt(viscosities[j]))) * pow(E, A);  
        }
    }
    return viscosityMixture = 1 / fluidity;
}
