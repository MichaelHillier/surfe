#include <modelling_methods_builder.h>
#include <modeling_methods.h>
#include <continuous_property.h>
#include <lajaunie.h>
#include <single_surface.h>
#include <stratigraphic_surfaces.h>
using namespace Surfe;

GRBF_Builder::GRBF_Builder() {}

GRBF_Modelling_Methods* GRBF_Builder::get_method(
    const model_parameters& m_parameters, const Constraints& input) {
    if (m_parameters.model_type == Parameter_Types::Single_surface)
        return new Single_Surface(m_parameters, input);
    else if (m_parameters.model_type == Parameter_Types::Lajaunie_approach)
        return new Lajaunie_Approach(m_parameters, input);
    else if (m_parameters.model_type == Parameter_Types::Stratigraphic_horizons)
        return new Stratigraphic_Surfaces(m_parameters, input);
    else
        return new Continuous_Property(m_parameters, input);
}
