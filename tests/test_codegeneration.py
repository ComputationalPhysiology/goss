import goss


def test_GossCodeGenerator(oscilator_ode):

    code_params = {
        "state_repr": "named",
        "body_repr": "named",
        "use_cse": False,
        "generate_forward_backward_subst": True,
        "generate_jacobian": True,
        "generate_lu_factorization": True,
        "optimize_exprs": "numerals",
    }
    codegen = goss.codegeneration.GossCodeGenerator(
        oscilator_ode,
        monitored=["energy"],
        field_states=["x"],
        field_parameters=["a"],
        code_params=code_params,
    )

    class_code = codegen.class_code()
    file_code = codegen.file_code()
    assert "#include <goss/ParameterizedODE.h>" in file_code
    assert "#include <goss/ParameterizedODE.h>" not in class_code
    assert class_code in file_code
    assert "void forward_backward_subst" in class_code
    assert "void lu_factorize" in class_code
    assert "void compute_jacobian" in class_code
    assert "jac[1] = -a;" in class_code
    assert "void linearized_eval" in class_code
    assert "monitored[0] = a*y + b*x;" in class_code
