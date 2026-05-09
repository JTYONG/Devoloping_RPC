#include<iostream>
#include<string>


int main(){
    const std::string path = std::getenv("GARFIELD_INSTALL");
    std::cout<<path<<std::endl;
    std::cout<<path+"/share/Garfield/Data/IonMobility_SF6-_SF6.txt"<<std::endl;
    return 0;
}