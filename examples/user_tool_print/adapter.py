"""AI-adapted entry point: preserve the submitted print inside run."""

def run(inputs, params, output_dir):
    print("ASTRAFEX_USER_PRINT_OK")
    return {"artifacts": [], "summary": {"science_validated": False}}
