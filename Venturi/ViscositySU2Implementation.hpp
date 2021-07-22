#include <iostream>
#include <cmath>
#include <math.h>
#include <vector>
#include <numeric>

class Viscosity {
    private:
        double viscosityMixture = 0.0;    
        std::vector<double> moleFractions; 
        std::vector<double> viscositySutherland; 
        double temperature; 
    public: 
        std::vector<double> massToMoleFractions(std::vector<double> massFractions, std::vector<double> &molarMasses, int numberOfComponents);
        std::vector<double> sutherlandViscosity(double temperature, std::vector<double> viscosityRef, std::vector<double> temperatureRef, std::vector<double> sutherlandCoeffRef, int numberOfComponents);
        double wilkeViscosity(std::vector<double> massFractions, std::vector<double> &molarMasses, std::vector<double> &viscositySutherland, int numberOfComponents);
        double davidsonViscosity(std::vector<double> massFractions, std::vector<double> &molarMasses, std::vector<double> &viscositySutherland, int numberOfComponents);
}; 

std::vector<double> Viscosity::massToMoleFractions(std::vector<double> massFractions, std::vector<double> &molarMasses, int numberOfComponents){
    double mixtureMolarMass = 0.0; 
    double mixtureFractionDenumerator = 0.0; 
    moleFractions.clear(); 

    for(int i = 0; i < numberOfComponents; i++){
        mixtureMolarMass += massFractions[i] / molarMasses[i]; 
    }

    for(int i = 0; i < numberOfComponents; i++){
        moleFractions.push_back((massFractions[i] / molarMasses[i]) / mixtureMolarMass);
        //std::cout << moleFractions[i] << " ";
    }
    //std::cout << "\n";
    return moleFractions; 
}

std::vector<double> Viscosity::sutherlandViscosity(double temperature, std::vector<double> viscosityRef, std::vector<double> temperatureRef, std::vector<double> sutherlandCoeffRef, int numberOfComponents){
    double temperatureNonDim = 0.0;  
    viscositySutherland.clear(); 

    for(int i = 0; i < numberOfComponents; i++){
        temperatureNonDim = temperature / temperatureRef[i];
        viscositySutherland.push_back(viscosityRef[i]*temperatureNonDim *sqrt(temperatureNonDim)*(temperatureRef[i]+sutherlandCoeffRef[i])/(temperature+sutherlandCoeffRef[i])); 
        std::cout << viscositySutherland[i] << " ";
    }
    std::cout << "\n";
    return viscositySutherland; 
}

double Viscosity::wilkeViscosity(std::vector<double> massFractions, std::vector<double> &molarMasses, std::vector<double> &viscositySutherland, int numberOfComponents){
    std::vector<double> phi; 
    std::vector<double> wilkeNumerator;
    std::vector<double> wilkeDenumeratorSum;
    double wilkeDenumerator = 0.0;
    viscosityMixture = 0.0; 
    wilkeDenumeratorSum.clear();
    wilkeNumerator.clear(); 
    
    massToMoleFractions(massFractions, molarMasses, numberOfComponents); 

    for(int i = 0; i < numberOfComponents; i++){
        for(int j = 0; j < numberOfComponents; j++){
            phi.push_back(pow((1 + pow(viscositySutherland[i] / viscositySutherland[j],0.5) * pow(molarMasses[j] / molarMasses[i],0.25)),2) / pow(8 * (1 + molarMasses[i] / molarMasses[j]),0.5)); 
            wilkeDenumerator += moleFractions[j] * phi[j];   
        }
        wilkeDenumeratorSum.push_back(wilkeDenumerator); 
        wilkeDenumerator = 0.0; 
        phi.clear();
        wilkeNumerator.push_back(moleFractions[i] * viscositySutherland[i]); 
        viscosityMixture += wilkeNumerator[i] / wilkeDenumeratorSum[i];  
    }
    return viscosityMixture; 
}

double Viscosity::davidsonViscosity(std::vector<double> massFractions, std::vector<double> &molarMasses, std::vector<double> &viscositySutherland, int numberOfComponents){
    double fluidity = 0.0;
    double E = 0.0; 
    double mixtureFractionDenumerator = 0.0;
    const double A = 0.375;   
    std::vector<double> mixtureFractions; 
    mixtureFractions.clear(); 
    
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
                fluidity += ((mixtureFractions[i] * mixtureFractions[j]) / (sqrt(viscositySutherland[i]) * sqrt(viscositySutherland[j]))) * pow(E, A);  
        }
    }
    return viscosityMixture = 1 / fluidity;
}

