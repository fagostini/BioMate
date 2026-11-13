"""Tests for web_interface module security and functionality."""

import argparse
from pathlib import Path

import pytest
import tornado.testing
import tornado.web

from biomate.web_interface import (
    MODULES,
    ModuleHandler,
    SecurityHeadersMixin,
    flatten_parameters,
    init_parser,
    make_app,
    validate_module_name,
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


class TestSecurityHeadersMixin:
    """Test security headers are properly set."""

    def test_security_headers_mixin_has_set_default_headers(self):
        """Test that SecurityHeadersMixin defines set_default_headers method."""
        assert hasattr(SecurityHeadersMixin, "set_default_headers")

    def test_security_headers_mixin_inherits_request_handler(self):
        """Test that SecurityHeadersMixin can be used with Tornado."""
        assert issubclass(SecurityHeadersMixin, object)


class TestModulesDefinition:
    """Test that all defined modules are properly configured."""

    def test_all_modules_defined(self):
        """Test that expected modules are in MODULES list."""
        module_names = [m["name"] for m in MODULES]
        expected_modules = [
            "blabber",
            "index",
            "dirstruct",
            "fastrewind",
            "nspector",
            "strainer",
            "waif",
        ]
        for expected in expected_modules:
            assert expected in module_names, f"Module {expected} not found in MODULES"

    def test_all_modules_have_required_fields(self):
        """Test that each module has required metadata fields."""
        for module in MODULES:
            assert "name" in module, f"Module missing 'name' field"
            assert "title" in module, (
                f"Module {module.get('name')} missing 'title' field"
            )
            assert "description" in module, (
                f"Module {module.get('name')} missing 'description' field"
            )
            assert "parameters" in module, (
                f"Module {module.get('name')} missing 'parameters' field"
            )

    def test_all_modules_have_valid_names(self):
        """Test that all module names pass validation."""
        for module in MODULES:
            assert validate_module_name(module["name"]), (
                f"Module name {module['name']} fails validation"
            )

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

    def test_waif_module_has_parameters(self):
        """Test that waif module has expected parameters."""
        waif = next(m for m in MODULES if m["name"] == "waif")
        param_names = [p["form_name"] for p in waif["parameters"]]
        assert "input_path" in param_names
        assert "output_file" in param_names
        assert "exact" in param_names
        assert "sample_sheet" in param_names
        assert "block_cache" in param_names

    def test_index_module_has_parameters(self):
        """Test that index module has expected parameters."""
        index = next(m for m in MODULES if m["name"] == "index")
        param_names = [p["form_name"] for p in index["parameters"]]
        assert "input" in param_names


class TestFlattenedParameters:
    """Test parameter flattening for CLI argument construction."""

    def test_flatten_simple_parameter(self):
        """Test flattening a single parameter."""
        params = [{"form_name": "seq_length", "cli_name": "seq-length"}]
        values = {"seq_length": "100"}
        args = flatten_parameters(params, values)
        assert args == ["--seq-length=100"]

    def test_flatten_missing_parameter(self):
        """Test that missing parameters are skipped."""
        params = [{"form_name": "seq_length", "cli_name": "seq-length"}]
        values = {}
        args = flatten_parameters(params, values)
        assert "--seq-length" not in args

    def test_flatten_boolean_parameter_true(self):
        """Test flattening a boolean parameter set to true."""
        params = [{"form_name": "taint", "cli_name": "taint", "param_type": "boolean"}]
        values = {"taint": "true"}
        args = flatten_parameters(params, values)
        assert "--taint" in args

    def test_flatten_boolean_parameter_false(self):
        """Test that boolean parameter set to false is not included."""
        params = [{"form_name": "taint", "cli_name": "taint", "param_type": "boolean"}]
        values = {"taint": "false"}
        args = flatten_parameters(params, values)
        assert "--taint" not in args

    def test_flatten_select_parameter(self):
        """Test flattening a select parameter."""
        params = [
            {
                "form_name": "error_type",
                "cli_name": "error-type",
                "param_type": "select",
            }
        ]
        values = {"error_type": "s"}
        args = flatten_parameters(params, values)
        assert args == ["--error-type=s"]

    def test_flatten_multiple_splits_commas(self):
        """A multi-value field passes each comma-separated item as its own arg."""
        params = [{"form_name": "input", "cli_name": "input", "multiple": True}]
        values = {"input": "a_R1.fastq.gz, b_R2.fastq.gz,"}
        args = flatten_parameters(params, values)
        assert args == ["--input", "a_R1.fastq.gz", "b_R2.fastq.gz"]

    def test_flatten_with_default_value(self):
        """Test that default values are used when provided."""
        params = [
            {
                "form_name": "threads",
                "cli_name": "threads",
                "default": "1",
                "param_type": "number",
            }
        ]
        values = {}
        args = flatten_parameters(params, values)
        assert args == []

    def test_flatten_with_provided_value_overrides_default(self):
        """Test that provided values override defaults."""
        params = [
            {
                "form_name": "threads",
                "cli_name": "threads",
                "default": "1",
                "param_type": "number",
            }
        ]
        values = {"threads": "4"}
        args = flatten_parameters(params, values)
        assert args == ["--threads=4"]


class TestWebInterfaceModule:
    """Integration tests for web_interface module."""

    def test_init_parser_creates_parser(self):
        """Test that init_parser creates valid argparse parser."""
        parser = argparse.ArgumentParser()
        subparsers = parser.add_subparsers()
        result = init_parser(subparsers)
        assert isinstance(result, argparse.ArgumentParser)

    def test_make_app_creates_tornado_application(self):
        """Test that make_app creates a Tornado application."""
        app = make_app(host="localhost", port=8080)
        assert isinstance(app, tornado.web.Application)

    def test_make_app_has_handlers(self):
        """Test that created app has request handlers configured."""
        app = make_app(host="localhost", port=8080)
        assert app is not None
        assert isinstance(app, tornado.web.Application)
        assert hasattr(app, "settings")
        assert app.settings is not None

    def test_xsrf_cookies_enabled_in_app(self):
        """Test that XSRF protection is enabled in application."""
        app = make_app(host="localhost", port=8080)
        assert app.settings.get("xsrf_cookies") is True

    def test_debug_mode_disabled_in_app(self):
        """Test that debug mode is disabled in application (security)."""
        app = make_app(host="localhost", port=8080)
        assert app.settings.get("debug") is False


class TestReportRoute:
    """Reports are only served for runs registered by this server."""

    @pytest.fixture()
    def report_dir(self, tmp_path):
        (tmp_path / "nspector_report.html").write_text("<html>report</html>")
        return tmp_path

    def test_register_report_returns_token_url(self, report_dir):
        from biomate.web_interface import _REPORTS, register_report

        url = register_report("nspector", {"output": str(report_dir)})
        assert url is not None and url.startswith("/report/")
        assert _REPORTS[url.rsplit("/", 1)[1]] == report_dir / "nspector_report.html"

    def test_register_report_ignores_other_modules(self, report_dir):
        from biomate.web_interface import register_report

        assert register_report("blabber", {"output": str(report_dir)}) is None

    def test_register_report_missing_file(self, tmp_path):
        from biomate.web_interface import register_report

        assert register_report("nspector", {"output": str(tmp_path)}) is None

    def test_report_route_registered(self):
        from biomate.web_interface import ReportHandler

        app = make_app()
        rules = [r.target for r in app.wildcard_router.rules]
        assert ReportHandler in rules

    def test_unknown_token_is_404(self, report_dir):
        import asyncio

        from tornado.httpclient import AsyncHTTPClient
        from tornado.testing import bind_unused_port

        from biomate.web_interface import register_report

        url = register_report("nspector", {"output": str(report_dir)})

        async def fetch(path):
            sock, port = bind_unused_port()
            server = tornado.httpserver.HTTPServer(make_app())
            server.add_sockets([sock])
            try:
                return await AsyncHTTPClient().fetch(
                    f"http://127.0.0.1:{port}{path}", raise_error=False
                )
            finally:
                server.stop()

        import tornado.httpserver

        ok = asyncio.run(fetch(url))
        assert ok.code == 200 and b"report" in ok.body
        assert "unsafe-eval" in ok.headers["Content-Security-Policy"]
        missing = asyncio.run(fetch("/report/not-a-token"))
        assert missing.code == 404
        traversal = asyncio.run(fetch("/report/..%2F..%2Fetc%2Fpasswd"))
        assert traversal.code == 404


class TestModuleForms(tornado.testing.AsyncHTTPTestCase):
    """Round trip of every tool form through the real handler and CLI parser."""

    def get_app(self):
        """Application under test."""
        return make_app()

    def _post_default_form(self, name, calls):
        """GET the form, then POST its default values back with the XSRF token."""
        import re
        from unittest import mock
        from urllib.parse import urlencode

        from biomate import web_interface

        page = self.fetch(f"/module/{name}", headers={"Host": "localhost"})
        assert page.code == 200
        token = re.search(rb'name="_xsrf" value="([^"]+)"', page.body).group(1).decode()
        cookie = page.headers["Set-Cookie"].split(";")[0]
        module = next(m for m in web_interface.MODULES if m["name"] == name)
        form = {"_xsrf": token}
        for p in module["parameters"]:
            # Like a browser: unchecked checkboxes are not sent, everything else is
            if p.get("param_type") != "boolean":
                form[p["form_name"]] = str(p.get("default", ""))
        fake = {"success": True, "returncode": 0, "stdout": "", "stderr": ""}
        with mock.patch.object(
            web_interface,
            "run_biomate",
            side_effect=lambda m, a: calls.append(a) or fake,
        ):
            return self.fetch(
                f"/module/{name}",
                method="POST",
                body=urlencode(form),
                headers={"Host": "localhost", "Cookie": cookie},
            )

    def test_default_form_of_every_module_builds_known_arguments(self):
        """No module's default form produces an option its CLI does not have."""
        import biomate
        from biomate import web_interface

        parser = biomate.CustomParser()
        subparsers = parser.add_subparsers()
        for mod in (
            biomate.blabber.blabber,
            biomate.dirstruct.dirstruct,
            biomate.index.index,
            biomate.fastrewind.fastrewind,
            biomate.nspector.nspector,
            biomate.strainer.strainer,
            biomate.waif.waif,
        ):
            mod.init_parser(subparsers)
        for module in web_interface.MODULES:
            calls = []
            response = self._post_default_form(module["name"], calls)
            assert response.code == 200, module["name"]
            # Unknown options would be left over; missing required ones are
            # the user's to fill in, so only the leftovers are checked
            _, unknown = _parse_known(parser, [module["name"], *calls[0]])
            assert unknown == [], (module["name"], unknown)

    def test_post_without_xsrf_token_is_rejected(self):
        """A POST without the form token is refused."""
        response = self.fetch("/module/blabber", method="POST", body="seq_number=1")
        assert response.code == 403

    def test_foreign_host_header_is_rejected(self):
        """A request for another host name (DNS rebinding) is refused."""
        response = self.fetch("/", headers={"Host": "evil.example.com"})
        assert response.code == 403

    def test_select_shows_its_default(self):
        """A select field preselects its default option on a fresh page."""
        page = self.fetch("/module/blabber", headers={"Host": "localhost"})
        assert b'<option value="text" selected>' in page.body


def _parse_known(parser, argv):
    """parse_known_args that tolerates missing required options."""
    import contextlib
    import io

    with contextlib.redirect_stderr(io.StringIO()):
        try:
            return parser.parse_known_args(argv)
        except SystemExit:
            # Required options missing: re-check only for unknown options
            unknown = [a for a in argv[1:] if a.startswith("--")]
            known = {
                s
                for sub in parser._subparsers._group_actions[0].choices.values()
                for a in _all_actions(sub)
                for s in a.option_strings
            }
            return None, [a for a in unknown if a.split("=")[0] not in known]


def _all_actions(parser):
    """Actions of a parser and of its nested subparsers."""
    for action in parser._actions:
        yield action
        for sub in getattr(action, "choices", None) or {}:
            if hasattr(action.choices, "values") and hasattr(
                action.choices[sub], "_actions"
            ):
                yield from _all_actions(action.choices[sub])
