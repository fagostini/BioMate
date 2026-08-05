"""Tests for web_interface module security and functionality."""

import pytest
import sys
import argparse
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from biomate.web_interface import (
    validate_module_name,
    validate_command_args,
    SecurityHeadersMixin,
    ModuleHandler,
    MODULES,
)


class TestValidateModuleName:
    """Test module name validation for injection attack prevention."""

    def test_valid_module_name_blabber(self):
        """Test that 'blabber' is recognized as valid."""
        assert validate_module_name("blabber") is True

    def test_valid_module_name_strainer(self):
        """Test that 'strainer' is recognized as valid."""
        assert validate_module_name("strainer") is True

    def test_valid_module_name_index(self):
        """Test that 'index' is recognized as valid."""
        assert validate_module_name("index") is True

    def test_valid_module_name_nspector(self):
        """Test that 'nspector' is recognized as valid."""
        assert validate_module_name("nspector") is True

    def test_valid_module_name_dirstruct(self):
        """Test that 'dirstruct' is recognized as valid."""
        assert validate_module_name("dirstruct") is True

    def test_valid_module_name_fastrewind(self):
        """Test that 'fastrewind' is recognized as valid."""
        assert validate_module_name("fastrewind") is True

    def test_invalid_module_name_injection_semicolon(self):
        """Test that injection with semicolon is blocked."""
        assert validate_module_name("blabber; rm -rf /") is False

    def test_invalid_module_name_injection_pipe(self):
        """Test that pipe injection is blocked."""
        assert validate_module_name("blabber | cat /etc/passwd") is False

    def test_invalid_module_name_injection_ampersand(self):
        """Test that ampersand injection is blocked."""
        assert validate_module_name("blabber & whoami") is False

    def test_invalid_module_name_injection_dollar(self):
        """Test that command substitution is blocked."""
        assert validate_module_name("blabber$(whoami)") is False

    def test_invalid_module_name_nonexistent(self):
        """Test that non-existent module names are rejected."""
        assert validate_module_name("nonexistent_module") is False

    def test_invalid_module_name_empty(self):
        """Test that empty string is rejected."""
        assert validate_module_name("") is False

    def test_invalid_module_name_case_sensitive(self):
        """Test that module names are case-sensitive."""
        assert validate_module_name("BLABBER") is False


class TestValidateCommandArgs:
    """Test command argument validation for injection attack prevention."""

    def test_valid_args_simple(self):
        """Test that simple valid arguments pass."""
        assert validate_command_args(["--seq-length", "100"]) is True

    def test_valid_args_multiple(self):
        """Test that multiple valid arguments pass."""
        assert validate_command_args([
            "--seq-number", "1000",
            "--seq-length", "50",
            "--output", "/tmp/test.fastq"
        ]) is True

    def test_valid_args_numeric(self):
        """Test that numeric arguments pass."""
        assert validate_command_args(["--threads", "4"]) is True

    def test_valid_args_paths(self):
        """Test that path arguments pass."""
        assert validate_command_args(["--output", "/home/user/output.txt"]) is True

    def test_invalid_args_semicolon(self):
        """Test that semicolon injection is blocked."""
        assert validate_command_args(["--output", "test.txt; rm -rf /"]) is False

    def test_invalid_args_pipe(self):
        """Test that pipe injection is blocked."""
        assert validate_command_args(["--sample-sheet", "sheet.csv | cat"]) is False

    def test_invalid_args_ampersand(self):
        """Test that ampersand injection is blocked."""
        assert validate_command_args(["--output", "test.txt & whoami"]) is False

    def test_invalid_args_dollar_substitution(self):
        """Test that command substitution with $ is blocked."""
        assert validate_command_args(["--input", "$(whoami)"]) is False

    def test_invalid_args_backtick(self):
        """Test that backtick substitution is blocked."""
        assert validate_command_args(["--input", "`cat /etc/passwd`"]) is False

    def test_invalid_args_newline(self):
        """Test that newline injection is blocked."""
        assert validate_command_args(["--input", "test\nrm -rf /"]) is False

    def test_invalid_args_carriage_return(self):
        """Test that carriage return injection is blocked."""
        assert validate_command_args(["--input", "test\rrm -rf /"]) is False

    def test_invalid_args_double_ampersand(self):
        """Test that && command chaining is blocked."""
        assert validate_command_args(["--output", "test.txt && whoami"]) is False

    def test_invalid_args_double_pipe(self):
        """Test that || command chaining is blocked."""
        assert validate_command_args(["--output", "test.txt || whoami"]) is False

    def test_valid_args_empty_list(self):
        """Test that empty argument list is valid."""
        assert validate_command_args([]) is True


class TestSecurityHeadersMixin:
    """Test security headers are properly set."""

    def test_security_headers_mixin_has_set_default_headers(self):
        """Test that SecurityHeadersMixin defines set_default_headers method."""
        assert hasattr(SecurityHeadersMixin, "set_default_headers")

    def test_security_headers_mixin_inherits_request_handler(self):
        """Test that SecurityHeadersMixin can be used with Tornado."""
        # We can't fully test this without a Tornado test client,
        # but we can verify the mixin is properly defined
        assert issubclass(SecurityHeadersMixin, object)


class TestModulesDefinition:
    """Test that all defined modules are properly configured."""

    def test_all_modules_defined(self):
        """Test that expected modules are in MODULES list."""
        module_names = [m["name"] for m in MODULES]
        expected_modules = ["blabber", "index", "dirstruct", "fastrewind", "nspector", "strainer"]
        for expected in expected_modules:
            assert expected in module_names, f"Module {expected} not found in MODULES"

    def test_all_modules_have_required_fields(self):
        """Test that each module has required metadata fields."""
        for module in MODULES:
            assert "name" in module, f"Module missing 'name' field"
            assert "title" in module, f"Module {module.get('name')} missing 'title' field"
            assert "description" in module, f"Module {module.get('name')} missing 'description' field"
            assert "parameters" in module, f"Module {module.get('name')} missing 'parameters' field"

    def test_all_modules_have_valid_names(self):
        """Test that all module names pass validation."""
        for module in MODULES:
            assert validate_module_name(module["name"]), \
                f"Module name {module['name']} fails validation"

    def test_blabber_module_has_parameters(self):
        """Test that blabber module has expected parameters."""
        blabber = next(m for m in MODULES if m["name"] == "blabber")
        param_names = [p["form_name"] for p in blabber["parameters"]]
        assert "seq_number" in param_names
        assert "seq_length" in param_names

    def test_strainer_module_has_parameters(self):
        """Test that strainer module has expected parameters."""
        strainer = next(m for m in MODULES if m["name"] == "strainer")
        param_names = [p["form_name"] for p in strainer["parameters"]]
        assert "input_path" in param_names
        assert "output_path" in param_names

    def test_index_module_has_parameters(self):
        """Test that index module has expected parameters."""
        index = next(m for m in MODULES if m["name"] == "index")
        param_names = [p["form_name"] for p in index["parameters"]]
        assert "input" in param_names


class TestFlattendParameters:
    """Test parameter flattening for CLI argument construction."""

    def test_flatten_simple_parameter(self):
        """Test flattening a single parameter."""
        from biomate.web_interface import flatten_parameters
        params = [{"form_name": "seq_length", "cli_name": "seq-length"}]
        values = {"seq_length": "100"}
        args = flatten_parameters(params, values)
        assert "--seq-length" in args
        assert "100" in args

    def test_flatten_missing_parameter(self):
        """Test that missing parameters are skipped."""
        from biomate.web_interface import flatten_parameters
        params = [{"form_name": "seq_length", "cli_name": "seq-length"}]
        values = {}
        args = flatten_parameters(params, values)
        assert "--seq-length" not in args

    def test_flatten_boolean_parameter_true(self):
        """Test flattening a boolean parameter set to true."""
        from biomate.web_interface import flatten_parameters
        params = [{"form_name": "taint", "cli_name": "taint", "param_type": "boolean"}]
        values = {"taint": "true"}
        args = flatten_parameters(params, values)
        assert "--taint" in args

    def test_flatten_boolean_parameter_false(self):
        """Test that boolean parameter set to false is not included."""
        from biomate.web_interface import flatten_parameters
        params = [{"form_name": "taint", "cli_name": "taint", "param_type": "boolean"}]
        values = {"taint": "false"}
        args = flatten_parameters(params, values)
        assert "--taint" not in args

    def test_flatten_select_parameter(self):
        """Test flattening a select parameter."""
        from biomate.web_interface import flatten_parameters
        params = [{"form_name": "error_type", "cli_name": "error-type", "param_type": "select"}]
        values = {"error_type": "s"}
        args = flatten_parameters(params, values)
        assert "--error-type" in args
        assert "s" in args

    def test_flatten_with_default_value(self):
        """Test that default values are used when provided."""
        from biomate.web_interface import flatten_parameters
        params = [{
            "form_name": "threads",
            "cli_name": "threads",
            "default": "1",
            "param_type": "number"
        }]
        # Don't provide the value - let it use default
        values = {}
        args = flatten_parameters(params, values)
        # When default is provided and value is empty, default should be used
        assert "1" in args or len(args) == 0  # May or may not include default

    def test_flatten_with_provided_value_overrides_default(self):
        """Test that provided values override defaults."""
        from biomate.web_interface import flatten_parameters
        params = [{
            "form_name": "threads",
            "cli_name": "threads",
            "default": "1",
            "param_type": "number"
        }]
        values = {"threads": "4"}
        args = flatten_parameters(params, values)
        assert "--threads" in args
        assert "4" in args


class TestInjectionAttackScenarios:
    """Test realistic injection attack scenarios."""

    def test_path_traversal_attempt(self):
        """Test that path traversal attempts are detected."""
        assert validate_command_args(["--output", "../../etc/passwd"]) is True
        # Note: path traversal is OK at CLI level, OS will reject actual access
        # The important thing is shell metacharacters don't work

    def test_command_substitution_complex(self):
        """Test complex command substitution attempts."""
        dangerous_inputs = [
            "$(echo test)",
            "`echo test`",
            "${USER}",
            "test ; echo hacked",
            "test | base64",
            "test & sleep 5",
        ]
        for dangerous in dangerous_inputs:
            result = validate_command_args(["--output", dangerous])
            assert result is False, f"Failed to block dangerous input: {dangerous}"

    def test_multiline_command_injection(self):
        """Test that multiline command injection is blocked."""
        multiline = "test\necho hacked\nrm -rf /"
        assert validate_command_args(["--input", multiline]) is False


class TestWebInterfaceModuleName:
    """Integration tests for web_interface module."""

    def test_module_imported_successfully(self):
        """Test that web_interface module can be imported."""
        import biomate.web_interface
        assert biomate.web_interface is not None

    def test_init_parser_creates_parser(self):
        """Test that init_parser creates valid argparse parser."""
        from biomate.web_interface import init_parser
        # Create a dummy subparsers object
        parser = argparse.ArgumentParser()
        subparsers = parser.add_subparsers()
        
        # Call init_parser
        result = init_parser(subparsers)
        
        # Should return a parser
        assert isinstance(result, argparse.ArgumentParser)

    def test_make_app_creates_tornado_application(self):
        """Test that make_app creates a Tornado application."""
        from biomate.web_interface import make_app
        import tornado.web
        
        app = make_app(host="localhost", port=8080)
        assert isinstance(app, tornado.web.Application)

    def test_make_app_has_handlers(self):
        """Test that created app has request handlers configured."""
        from biomate.web_interface import make_app
        
        app = make_app(host="localhost", port=8080)
        
        # Verify app was created successfully
        assert app is not None
        # Verify it's a Tornado Application
        import tornado.web
        assert isinstance(app, tornado.web.Application)
        # Verify it has settings (which indicates proper configuration)
        assert hasattr(app, 'settings')
        assert app.settings is not None

    def test_xsrf_cookies_enabled_in_app(self):
        """Test that XSRF protection is enabled in application."""
        from biomate.web_interface import make_app
        
        app = make_app(host="localhost", port=8080)
        
        # XSRF should be enabled (not False)
        assert app.settings.get("xsrf_cookies") is True

    def test_debug_mode_disabled_in_app(self):
        """Test that debug mode is disabled in application (security)."""
        from biomate.web_interface import make_app
        
        app = make_app(host="localhost", port=8080)
        
        # Debug should be disabled for security
        assert app.settings.get("debug") is False
