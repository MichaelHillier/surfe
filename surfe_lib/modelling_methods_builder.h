#ifndef MODELLING_METHODS_BUILDER_H
#define MODELLING_METHODS_BUILDER_H

#include <modelling_input.h>
#include <modelling_parameters.h>
#include <modeling_methods.h>

class SURFE_LIB_EXPORT GRBF_Builder {
public:
    GRBF_Builder();
    GRBF_Modelling_Methods* get_method(const model_parameters& m_paramers,
                                       const Constraints& input);
};

#endif /* MODELLING_METHODS_BUILDER_H */
