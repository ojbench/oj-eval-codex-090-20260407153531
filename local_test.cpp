#include <iostream>
#define STRIFY(x) #x
#include STRIFY(fraction.hpp)
#include STRIFY(src.hpp)

const int interface_size = 3; const int connection_size = 3;
int from_[connection_size] = {1,1,2};
int to_[connection_size]   = {2,3,3};
fraction resistance_[connection_size] = {fraction(1,2), fraction(1,4), fraction(2)};
fraction current_[interface_size] = {fraction(2), fraction(1), fraction(-3)};
fraction voltage_[interface_size] = {fraction(1), fraction(2), fraction(1,2)};

int main(){
    resistive_network net(interface_size, connection_size, from_, to_, resistance_);
    std::cout << net.get_equivalent_resistance(1,2) << std::endl;
    std::cout << net.get_voltage(2, current_) << std::endl;
    std::cout << net.get_power(voltage_) << std::endl;
}
