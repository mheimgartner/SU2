#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iterator>
#include <algorithm>
#include "ViscosityGrid.hpp"

/* Input info
double W_ch4 = 16.0430;
double W_air = 28.9646; 
double visc_air = 1.845 * pow(10, -5);
double visc_ch4 = 1.113 * pow(10, -5);
*/

int main(){
    Viscosity Fluid; 
    std::vector< std::vector<double> > data;
    std::vector<double> viscosities;
    std::vector<double> molarMasses;
    std::vector<double> massFractions; 
    std::vector<double> YCH4; 
    std::vector<double> YAIR; 
    int numberOfComponents; 
    double viscosityValue;  
    double molarMassValue; 
    double value; 
    double massFractionMethane; 
    double massFractionAir; 

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

    std::ifstream MyReadFile("massFractionsLinspace.txt"); // Read .txt file containing mass fractions CH4 at outlet
    std::ofstream MyWriteFile("viscositiesMixture.csv"); // Write .csv file containing Wilke and Davidson viscosities for all grid points at the outlet
    
    if (MyReadFile.is_open()) { // If file has correctly opened...
		// Dynamically store data into 2D vector 
        while (MyReadFile.good()) { // ... and while there are no errors
            for(int row = 0; row < 1001; row++){ // Remember to change this value when the text file number of gridpoints is changed
                std::vector<double> temp; 
                MyReadFile >> massFractionMethane >> massFractionAir; 
                YCH4.push_back(massFractionMethane);
                YAIR.push_back(massFractionAir);
                for(int column = 0; column < numberOfComponents; column++){
                    temp.push_back(massFractionMethane);
                    temp.push_back(massFractionAir); 
                }
                data.push_back(temp);
            }
		}
    }
    else std::cout << "Unable to open file" << std::endl;

    MyWriteFile << "YCH4" << "," << "viscosityWilke" << "," 
    << "viscosityDavidson" << std::endl; 

    for(int i = 0; i < data.size(); i++){   // Loop over rows  
        for(int j = 0; j < data[i].size(); j++){ // Loop over columns
            massFractions.push_back(data[i][j]);
            //std::cout << massFractions[j] << " ";  
        }
        MyWriteFile << YCH4[i] << "," << Fluid.wilkeViscosity(massFractions, molarMasses, viscosities, numberOfComponents) << "," 
        << Fluid.davidsonViscosity(massFractions, molarMasses, viscosities, numberOfComponents) << std::endl; 
        massFractions.clear(); 
    } 

    //std::cout << std::endl;

    MyReadFile.close(); // Close the file

    return 0;
}
