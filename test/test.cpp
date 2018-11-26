#include <test.h>
#include <modelling_input.h>
#include <continuous_property.h>
#include <modeling_methods.h>
#include <vector>
using namespace Surfe;

int main(int argc, char* argv[]) {
    Basic_input input = Basic_input();
    std::vector<Interface> interfaces;
    std::vector<Evaluation_Point> eval;
    double z1 = 0.0;
    double z2 = 1.0;
    interfaces.push_back(Interface(0.0, 0.0, 0.0, 0.0));
    for (double x = 0.0; x < 10.00; x += 1) {
        for (double y = 0.0; y < 10.00; y += 1) {
            for (double z = 0.0; z < 5.00; z += 1) {
                eval.push_back(Evaluation_Point(x, y, z));
            }
        }
    }
    eval.push_back(Evaluation_Point(0.1, 1.1, 1.1));
    eval.push_back(Evaluation_Point(1.1, 1.1, 1.1));
    interfaces.push_back(Interface(1.0, 1.0, 1.0, 1.0));
    for (double x = 0.0; x < 10.0; x += 1.0) {
        for (double y = 0.0; y < 10.0; y += 1.0) {
            interfaces.push_back(Interface(x, y, z1, z1));
            interfaces.push_back(Interface(x, y, z2, z2));
        }
    }
    input.itrface = &interfaces;
    input.evaluation_pts = &eval;
    model_parameters parameters = model_parameters();
    Continuous_Property* metho = new Continuous_Property(parameters, input);
    metho->run_algorithm();
    std::cout << "Finished" << std::endl;
    delete metho;
}
