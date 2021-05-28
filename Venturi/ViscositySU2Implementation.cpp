#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "ViscositySU2Implementation.hpp"
//#include "../SU2_CFD/include/fluid/CSutherland.hpp"

int main(){
    Viscosity Fluid; 
    std::vector< std::vector<double> > data;
    //std::vector<double> viscositySutherland;
    std::vector<double> molarMasses;
    std::vector<double> massFractions; 
    std::vector<double> viscosityRef; 
    std::vector<double> temperatureRef; 
    std::vector<double> sutherlandCoeffRef; 
    std::vector<double> temperatureData;  
    std::vector<double> YCH4; 
    std::vector<double> YAIR; 
    int numberOfComponents; 
    double viscosityValue;  
    double molarMassValue; 
    double temperatureRefValue; 
    double viscosityRefValue; 
    double sRefValue; 
    double value; 
    double massFractionMethane; 
    double massFractionAir; 
    double temperatureValue;
    double temperature; 
    std::vector<double> viscositySutherland; 

    std::cout << "Enter number of components in mixture: "; 
    std::cin >> numberOfComponents; 
    
    for(int i = 0; i < numberOfComponents; i++){
        std::cout << "Enter mu_ref, t_ref, s_ref and the molar mass of component " << i+1 <<  ":" << std::endl;
        std::cin >> viscosityRefValue;
        viscosityRef.push_back(viscosityRefValue);

        std::cin >> temperatureRefValue;
        temperatureRef.push_back(temperatureRefValue);

        std::cin >> sRefValue;
        sutherlandCoeffRef.push_back(sRefValue);

        std::cin >> molarMassValue;
        molarMasses.push_back(molarMassValue); 
    }

    std::cout << std::endl;

    for(int i = 0; i < numberOfComponents; i++){
        std::cout << "Component " << i+1 <<  ":" << std::endl;
        std::cout << "mu_ref: " << viscosityRef[i] << ", t_ref: " << temperatureRef[i] << ", s_ref: " 
        << sutherlandCoeffRef[i] << ", Molar mass: " << molarMasses[i] << std::endl;
        std::cout << std::endl;
    }

    std::ifstream MyReadFile("inputDataT1100.txt"); // Read .txt file containing mass fractions CH4 at outlet
    std::ofstream MyWriteFile("viscositiesMixtureBoushehriExactFitFullSutherlandT1100.csv"); // Write .csv file containing Wilke and Davidson viscosities for all grid points at the outlet
    
    if (MyReadFile.is_open()) { // If file has correctly opened...
		// Dynamically store data into 2D vector 
        while (MyReadFile.good()) { // ... and while there are no errors
            for(int row = 0; row < 1001; row++){ // Remember to change this value when the text file number of gridpoints is changed
                std::vector<double> temporaryMassFractions; 
                MyReadFile >> temperatureValue >> massFractionMethane >> massFractionAir; 
                temperatureData.push_back(temperatureValue); 
                YCH4.push_back(massFractionMethane);
                YAIR.push_back(massFractionAir);
                for(int column = 0; column < numberOfComponents; column++){
                    temporaryMassFractions.push_back(massFractionMethane);
                    temporaryMassFractions.push_back(massFractionAir); 
                }
                data.push_back(temporaryMassFractions);
            }
		}
    }
    else std::cout << "Unable to open file" << std::endl;

    //zorgen dat viscosities vanuit sutherland functie als input gebruikt worden in de viscositymixture functies 

    MyWriteFile << "Temperature" << "," << "YCH4" << "," << "viscosityWilke" << "," 
    << "viscosityDavidson" << std::endl; 

    for(int i = 0; i < data.size(); i++){   // Loop over rows  
        viscositySutherland = Fluid.sutherlandViscosity(temperatureData[i], viscosityRef, temperatureRef, sutherlandCoeffRef, numberOfComponents); 

        for(int j = 0; j < data[i].size(); j++){ // Loop over columns
            massFractions.push_back(data[i][j]);
            //std::cout << massFractions[j] << " "; 
            //std::cout << viscositySutherland[j] << " ";  
        }
        
        MyWriteFile << temperatureData[i] << "," << YCH4[i] << "," 
        << Fluid.wilkeViscosity(massFractions, molarMasses, viscositySutherland, numberOfComponents) << "," 
        << Fluid.davidsonViscosity(massFractions, molarMasses, viscositySutherland, numberOfComponents) << std::endl; 
      
        massFractions.clear(); 
    } 

    MyReadFile.close(); // Close the file

    return 0;
}
