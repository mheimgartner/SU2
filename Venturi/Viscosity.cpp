#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "Viscosity.hpp"

/* Input info
double W_ch4 = 16.0430;
double W_air = 28.9646; 
double visc_air = 1.845 * pow(10, -5);
double visc_ch4 = 1.113 * pow(10, -5);
*/

int main(){
    Viscosity Fluid; 
    std::vector<double> viscosities;
    std::vector<double> molarMasses;
    std::vector<double> massFractions; 
    int numberOfComponents; 
    double viscosityValue;  
    double molarMassValue; 
    double massFractionValue; 
    std::string massFractionOutletData; 

    std::cout << "Enter number of components in mixture: "; 
    std::cin >> numberOfComponents; 
    
    for(int i = 0; i < numberOfComponents; i++){
        std::cout << "Enter viscosity and molar mass of component " << i+1 <<  ":" << std::endl;
        std::cin >> viscosityValue;
        viscosities.push_back(viscosityValue);
        std::cin >> molarMassValue;
        molarMasses.push_back(molarMassValue); 
    }

    std::cout << std::endl;

    for(int i = 0; i < numberOfComponents; i++){
        std::cout << "Component " << i+1 <<  ":" << std::endl;
        std::cout << "Viscosity: " << viscosities[i] << ", Molar mass: " << molarMasses[i] << std::endl;
        std::cout << std::endl;
    }

    std::cout << std::endl;

    for(int i = 0; i < numberOfComponents; i++){
        std::cout << "Enter mass fraction of component " << i+1 <<  ":" << std::endl;
        std::cin >> massFractionValue;
        massFractions.push_back(massFractionValue);
    }

    std::cout << "Mixture viscosity Wilke: " << Fluid.wilkeViscosity(massFractions, molarMasses, viscosities, numberOfComponents) << " kg/(m s)" << std::endl;
    std::cout << "Mixture viscosity Davidson: " << Fluid.davidsonViscosity(massFractions, molarMasses, viscosities, numberOfComponents) << " kg/(m s)" << std::endl;
    
    std::ifstream MyReadFile("y_ch4.txt"); // Read .txt file containing mass fractions CH4 at outlet 

    // Use a while loop together with the getline() function to read the file line by line
    while (std::getline (MyReadFile, massFractionOutletData)){
        std::cout << massFractionOutletData << std::endl; 
    }

    MyReadFile.close(); // Close the file

    return 0;
}
