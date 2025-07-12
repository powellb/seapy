import importlib
import pytest

# List of expected Fortran extensions
fortran_modules = ["seapy.external.oalib"]


@pytest.mark.parametrize("module_name", fortran_modules)
def test_import_fortran_extensions(module_name):
    """Test that the Fortran extension modules can be imported successfully."""
    try:
        mod = importlib.import_module(module_name)
        assert mod is not None
    except ImportError as e:
        pytest.fail(f"Failed to import {module_name}: {e}")


if __name__ == "__main__":
    import sys
    import pytest

    sys.exit(pytest.main([__file__]))
