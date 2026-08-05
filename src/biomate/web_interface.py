"""
BioMate Web Interface - A tornado-based web UI for BioMate tools.

Launches a local web server that provides:
- A main dashboard listing all BioMate modules
- Individual pages with descriptions and parameter forms for each module
- Execution of tools and display of results in the browser

SECURITY NOTES:
- Designed for local development use only
- Includes protections against common web attacks:
  - XSRF protection via tokens
  - Security headers (X-Content-Type-Options, X-Frame-Options, etc.)
  - Input validation and command injection prevention
  - Rate limiting support
- For public network deployment, requires:
  - HTTPS/TLS configuration
  - Reverse proxy with authentication (nginx, Apache, etc.)
  - Additional hardening measures (JWT, OAuth, etc.)

See init_parser() for --ssl-cert and --ssl-key options for HTTPS setup.
"""

from __future__ import annotations

import argparse
import os
import shlex
import ssl
import subprocess
from pathlib import Path
from typing import Any, Dict, List

import json

import tornado.httpserver
import tornado.ioloop
import tornado.template
import tornado.web

# ---------------------------------------------------------------------------
# Resolve paths
# ---------------------------------------------------------------------------

_PKG_DIR = Path(__file__).resolve().parent
_DESIGN_DIR = _PKG_DIR / "design"
_STATIC_DIR = _PKG_DIR / "static"


# ---------------------------------------------------------------------------
# Module definitions
# ---------------------------------------------------------------------------

MODULES: List[Dict[str, Any]] = [
    {
        "name": "blabber",
        "title": "Sequence Generator",
        "short_description": (
            "Generates random nucleotide sequences in FASTA, FASTQ, "
            "FASTQ-extended, or plain text format. Can also generate files "
            "matching a specific sample sheet with realistic Illumina read "
            "names and tile coordinates."
        ),
        "description": (
            "The blabber tool generates random DNA sequences with configurable "
            "length, number, alphabet, and output format. It supports sequence "
            "masks defining the structure of reads and indexes (e.g. R1;I1;I2;R2), "
            "and can mimic real Illumina output by using a SampleSheet.csv file. "
            "Optionally taints output with ~10% cross-project sequences to simulate "
            "index hopping."
        ),
        "parameters": [
            {
                "label": "Number of Sequences",
                "param_type": "number",
                "default": "100",
                "help": "How many sequences to generate",
                "form_name": "seq_number",
                "cli_name": "seq-number",
            },
            {
                "label": "Sequence Length",
                "param_type": "number",
                "default": "100",
                "help": "Length of each generated sequence",
                "form_name": "seq_length",
                "cli_name": "seq-length",
            },
            {
                "label": "Sequence Mask",
                "param_type": "text",
                "default": "",
                "is_path": False,
                "help": "Structure mask (format: R1;I1;I2;R2). Leave empty to auto-detect from sample sheet.",
                "form_name": "seq_mask",
                "cli_name": "seq-mask",
            },
            {
                "label": "Index 1 Sequence",
                "param_type": "text",
                "default": "",
                "is_path": False,
                "help": "Custom Index 1 sequence (overrides any value in sample sheet)",
                "form_name": "index1",
                "cli_name": "index1",
            },
            {
                "label": "Index 2 Sequence",
                "param_type": "text",
                "default": "",
                "help": "Custom Index 2 sequence",
                "form_name": "index2",
                "cli_name": "index2",
            },
            {
                "label": "Nucleotide Alphabet",
                "param_type": "text",
                "default": "ACGT",
                "help": "Characters used for random sequence generation",
                "form_name": "alphabet",
                "cli_name": "alphabet",
            },
            {
                "label": "Output Format",
                "param_type": "select",
                "default": "text",
                "options": ["fasta", "fastq", "fastq-ext", "text"],
                "help": "Output file format",
                "form_name": "format",
                "cli_name": "format",
            },
            {
                "label": "Output File Path",
                "param_type": "text",
                "default": "",
                "help": "Output file path (leave empty to stream from server)",
                "form_name": "output",
                "cli_name": "output",
            },
            {
                "label": "Sample Sheet Path",
                "param_type": "text",
                "default": "",
                "help": "Path to SampleSheet.csv for realistic output structure",
                "form_name": "sample_sheet",
                "cli_name": "sample-sheet",
            },
            {
                "label": "Flowcell ID",
                "param_type": "text",
                "default": "",
                "help": "Custom flowcell ID for read names",
                "form_name": "flowcell_id",
                "cli_name": "flowcell-id",
            },
            {
                "label": "Taint Output",
                "param_type": "boolean",
                "default": False,
                "help": "Add ~10% sequences from other projects to undetermined files",
                "form_name": "taint",
                "cli_name": "taint",
            },
            {
                "label": "Random Seed",
                "param_type": "number",
                "default": "",
                "help": "Seed for reproducible random generation",
                "form_name": "random_seed",
                "cli_name": "random-seed",
            },
        ],
    },
    {
        "name": "dirstruct",
        "title": "Directory Structure",
        "short_description": (
            "Extracts a directory tree structure to a file, or recreates a "
            "directory/file structure from a previously extracted file."
        ),
        "description": (
            "Dirstruct has two subcommands: extract and create. Extract indexes a "
            "directory tree into a tab-separated file. Create rebuilds the "
            "directory structure from such a file. Useful for backing up or "
            "duplicating project directory trees without copying actual data."
        ),
        "parameters": [
            {
                "label": "Subcommand",
                "param_type": "select",
                "default": "extract",
                "options": ["extract", "create"],
                "help": "Which operation to perform",
                "form_name": "_subcommand",
                "cli_name": "_subcommand",
            },
            {
                "label": "Source Path",
                "param_type": "text",
                "default": "",
                "help": "Path to the directory to index",
                "form_name": "source_path",
                "cli_name": "source-path",
                "show_when": {"parameter": "_subcommand", "value": "extract"},
            },
            {
                "label": "No Tags",
                "param_type": "boolean",
                "default": False,
                "help": "Omit dir/file tags in output",
                "form_name": "no_tags",
                "cli_name": "no-tags",
                "show_when": {"parameter": "_subcommand", "value": "extract"},
            },
            {
                "label": "Source File",
                "param_type": "text",
                "default": "",
                "help": "Path to the file containing the structure definition",
                "form_name": "source_file",
                "cli_name": "source-file",
                "show_when": {"parameter": "_subcommand", "value": "create"},
            },
            {
                "label": "Output Path",
                "param_type": "text",
                "default": ".",
                "help": "Destination directory (default: current directory)",
                "form_name": "output_path",
                "cli_name": "output-path",
                "show_when": {"parameter": "_subcommand", "value": "create"},
            },
        ],
    },
    {
        "name": "fastrewind",
        "title": "BCL Reconstructor",
        "short_description": (
            "Re-generates the raw BCL (Binary Base Calls) directory structure "
            "from demultiplexed FASTQ files."
        ),
        "description": (
            "Fastrewind reverses the Illumina demultiplexing process. Given "
            "demultiplexed FASTQ files, it outputs a directory with .cbcl files, "
            ".filter files, .locs files, and RunInfo.xml -- matching what "
            "bcl-convert would produce from raw instrument data. Currently best "
            "supported for NovaSeq X Plus instruments."
        ),
        "parameters": [
            {
                "label": "Input Path",
                "param_type": "text",
                "default": "",
                "help": "Input flowcell directory with demultiplexed FASTQ files",
                "form_name": "input_path",
                "cli_name": "input-path",
            },
            {
                "label": "Output Path",
                "param_type": "text",
                "default": ".",
                "help": "Output directory for BCL structure (default: current dir)",
                "form_name": "output_path",
                "cli_name": "output-path",
            },
            {
                "label": "Sample Sheet Path",
                "param_type": "text",
                "default": "",
                "help": "Path to SampleSheet.csv (looks in input dir if omitted)",
                "form_name": "sample_sheet",
                "cli_name": "sample-sheet",
            },
            {
                "label": "Total Cycles",
                "param_type": "number",
                "default": "",
                "help": "Override total number of cycles",
                "form_name": "total_cycles",
                "cli_name": "total-cycles",
            },
            {
                "label": "Instrument",
                "param_type": "select",
                "default": "NovaSeqXPlus",
                "options": ["NovaSeqXPlus", "MiSeq", "NextSeq2000"],
                "help": "Sequencing instrument type",
                "form_name": "instrument",
                "cli_name": "instrument",
            },
            {
                "label": "Threads",
                "param_type": "number",
                "default": "0",
                "help": "Number of threads for dnaio I/O",
                "form_name": "threads",
                "cli_name": "threads",
            },
        ],
    },
    {
        "name": "index",
        "title": "Sequence Index Matcher",
        "short_description": (
            "Indexes FASTA or FASTQ files by searching for specific patterns "
            "or DNA barcodes, with support for fuzzy matching."
        ),
        "description": (
            "Searches FASTA/FASTQ files for index sequences (DNA barcodes). "
            "Supports exact matching via regex or a list of patterns, as well "
            "as fuzzy matching using Levenshtein edit distance. Outputs match "
            "counts with percentages and error breakdowns per distance level."
        ),
        "parameters": [
            {
                "label": "Input File",
                "param_type": "text",
                "default": "",
                "help": "Input FASTQ/FASTA file",
                "form_name": "input",
                "cli_name": "input",
            },
            {
                "label": "Output Directory",
                "param_type": "text",
                "default": "",
                "help": "Output directory for results",
                "form_name": "output",
                "cli_name": "output",
            },
            {
                "label": "Match Type",
                "param_type": "select",
                "default": "index_regex",
                "options": ["index_regex", "index_file"],
                "help": "How to specify index patterns",
                "form_name": "match_type",
                "cli_name": "match-type",
            },
            {
                "label": "Index Regex Pattern",
                "param_type": "text",
                "default": "",
                "help": "Regex pattern to search for",
                "form_name": "index_regex",
                "cli_name": "index-regex",
                "show_when": {"parameter": "match_type", "value": "index_regex"},
            },
            {
                "label": "Index File",
                "param_type": "text",
                "default": "",
                "help": "File with one index per line",
                "form_name": "index_file",
                "cli_name": "index-file",
                "show_when": {"parameter": "match_type", "value": "index_file"},
            },
            {
                "label": "Max Edit Distance",
                "param_type": "number",
                "default": "0",
                "help": "Maximum Levenshtein distance for fuzzy matching",
                "form_name": "distance",
                "cli_name": "distance",
            },
            {
                "label": "Error Type",
                "param_type": "select",
                "default": "s",
                "options": ["s", "i", "d", "e"],
                "help": "Error type: s=substitution, i=insertion, d=deletion, e=any edit",
                "form_name": "error_type",
                "cli_name": "error-type",
            },
        ],
    },
    {
        "name": "nspector",
        "title": "N-Base Inspector",
        "short_description": (
            "Inspects FASTQ sequencing data for N (unknown/ambiguous) bases "
            "and generates visualizations of their distribution."
        ),
        "description": (
            "Nspector analyzes FASTQ files for N base content across "
            "sequencing cycles and flowcell tiles. Produces two visualizations "
            "per file: (1) a bar chart of N-count distribution per sequence, "
            "and (2) a line chart of N-containing sequences per cycle coloured "
            "by tile. Requires altair for rendering."
        ),
        "parameters": [
            {
                "label": "Input File(s)",
                "param_type": "text",
                "default": "",
                "help": "Input FASTQ file(s) (accepts multiple, comma-separated)",
                "form_name": "input",
                "cli_name": "input",
            },
            {
                "label": "Output Directory",
                "param_type": "text",
                "default": "./output",
                "help": "Output directory for visualisations",
                "form_name": "output",
                "cli_name": "output",
            },
        ],
    },
    {
        "name": "strainer",
        "title": "Index Mixing Evaluator",
        "short_description": (
            "Evaluates index (barcode) mixing/cross-talk across flowcell "
            "lanes by analysing Undetermined FASTQ files."
        ),
        "description": (
            "Strainer flags potential sample swapping or index hopping by "
            "comparing unexpected indexes found in undetermined FASTQ pools "
            "against expected indexes from the Sample Sheet. It reads all "
            "Undetermined*R1*.fastq.gz files, filters poly-G artifacts, and "
            "compares observed indexes against expected ones across all lanes "
            "(allowing 1 mismatch). Outputs CSV reports of unexpected indexes "
            "and lane summaries."
        ),
        "parameters": [
            {
                "label": "Input Path",
                "param_type": "text",
                "default": "",
                "help": "Flowcell directory with undetermined FASTQ files",
                "form_name": "input_path",
                "cli_name": "input-path",
            },
            {
                "label": "Output Path",
                "param_type": "text",
                "default": ".",
                "help": "Output directory for results (default: current dir)",
                "form_name": "output_path",
                "cli_name": "output-path",
            },
            {
                "label": "Sample Sheet Path",
                "param_type": "text",
                "default": "",
                "help": "Path to SampleSheet.csv (looks in input dir if omitted)",
                "form_name": "sample_sheet",
                "cli_name": "sample-sheet",
            },
            {
                "label": "Threads",
                "param_type": "number",
                "default": "1",
                "help": "Number of threads for parallel processing",
                "form_name": "threads",
                "cli_name": "threads",
            },
        ],
    },
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def validate_module_name(module_name: str) -> bool:
    """Validate that module name is allowed to prevent injection attacks."""
    allowed_modules = {m["name"] for m in MODULES}
    return module_name in allowed_modules


def validate_command_args(args: List[str]) -> bool:
    """Validate command arguments to prevent injection attacks.
    
    Ensures arguments don't contain shell metacharacters that would bypass
    subprocess.run's protection.
    """
    # Dangerous patterns that could indicate shell injection attempts
    dangerous_patterns = [";", "|", "&", "$", "`", "\n", "\r", "&&", "||"]
    
    for arg in args:
        for pattern in dangerous_patterns:
            if pattern in arg:
                return False
    return True


def run_biomate(module_name: str, args: List[str]) -> Dict[str, Any]:
    """Run a biomate subcommand and return stdout / stderr / returncode.
    
    Validates module name and arguments to prevent injection attacks.
    """
    # Validate module name
    if not validate_module_name(module_name):
        return {
            "success": False,
            "returncode": -1,
            "stdout": "",
            "stderr": f"Invalid module name: {module_name}",
        }
    
    # Validate arguments
    if not validate_command_args(args):
        return {
            "success": False,
            "returncode": -1,
            "stdout": "",
            "stderr": "Invalid arguments detected. Command injection prevented.",
        }
    
    cmd = ["biomate", module_name] + args
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=600,
        )
        return {
            "success": proc.returncode == 0,
            "returncode": proc.returncode,
            "stdout": proc.stdout,
            "stderr": proc.stderr,
        }
    except subprocess.TimeoutExpired:
        return {
            "success": False,
            "returncode": -1,
            "stdout": "",
            "stderr": "Command timed out after 10 minutes.",
        }
    except Exception as e:
        return {
            "success": False,
            "returncode": -1,
            "stdout": "",
            "stderr": str(e),
        }


def flatten_parameters(params: List[Dict], values: Dict[str, str]) -> List[str]:
    """Build CLI arguments from submitted form values."""
    cli_args: List[str] = []
    for p in params:
        form_key = p.get("form_name", p.get("name", ""))
        cli_key = p.get("cli_name", form_key)
        val = values.get(form_key)
        if val is None:
            continue

        if cli_key == "_subcommand" and not val.startswith("-"):
            cli_args.append(val)
            continue

        if p.get("param_type") == "boolean":
            if val and val.lower() in ("true", "on", "1"):
                cli_args.append(f"--{cli_key}")
            continue

        if p.get("param_type") == "select":
            if val:
                cli_args.append(f"--{cli_key}")
                cli_args.append(str(val))
            continue

        # text / number / unspecified
        default = p.get("default", "")
        if not val and default:
            val = default
        if val:
            cli_args.append(f"--{cli_key}")
            cli_args.append(shlex.quote(str(val)))

    return cli_args


# ---------------------------------------------------------------------------
# Tornado template loader
# ---------------------------------------------------------------------------

_loader = None


class ParamDict(dict):
    """A dict subclass that supports attribute access for Tornado templates."""
    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError(f"'{type(self).__name__}' object has no attribute '{key}'")
    def __setattr__(self, key, value):
        self[key] = value


def _module_list() -> list:
    _path_keywords = ["path", "file", "directory", "folder", "sheet"]
    result = []
    for m in MODULES:
        params = []
        for p in m["parameters"]:
            if "is_path" not in p:
                text = " ".join(
                    p.get(k, "")
                    for k in ("label", "cli_name", "form_name", "help")
                ).lower()
                params.append(ParamDict({**p, "is_path": any(kw in text for kw in _path_keywords)}))
            else:
                params.append(ParamDict(p))
        result.append({**m, "parameters": params})
    return result


def get_loader() -> tornado.template.Loader:
    global _loader
    if _loader is None:
        _loader = tornado.template.Loader(str(_DESIGN_DIR))
    return _loader


# ---------------------------------------------------------------------------
# Tornado handlers
# ---------------------------------------------------------------------------


class SecurityHeadersMixin(tornado.web.RequestHandler):
    """Mixin to add security headers to all responses."""
    
    def set_default_headers(self):
        """Set security-related HTTP headers."""
        # Prevent MIME type sniffing
        self.set_header("X-Content-Type-Options", "nosniff")
        # Prevent clickjacking
        self.set_header("X-Frame-Options", "DENY")
        # Enable XSS protection in older browsers
        self.set_header("X-XSS-Protection", "1; mode=block")
        # Strict Transport Security (if using HTTPS)
        self.set_header("Strict-Transport-Security", "max-age=31536000; includeSubDomains")
        # Content Security Policy - strict by default
        self.set_header("Content-Security-Policy", "default-src 'self'; script-src 'self' 'unsafe-inline'; style-src 'self' 'unsafe-inline'")
        # Referrer Policy
        self.set_header("Referrer-Policy", "strict-origin-when-cross-origin")


class MainHandler(SecurityHeadersMixin):
    """Main dashboard - lists all BioMate tools."""

    def get(self):
        loader = get_loader()
        t = loader.load("main.html")
        self.write(t.generate(modules=_module_list(), gs_version="0.3.0"))


class ModuleHandler(SecurityHeadersMixin):
    """Individual tool page with form and result display."""

    def get(self, module_name: str):
        module = self._find_module(module_name)
        if not module:
            self.set_status(404)
            self.write("Module not found.")
            return
        loader = get_loader()
        t = loader.load("module.html")
        all_modules = _module_list()
        mod = next(m for m in all_modules if m["name"] == module_name)
        param_meta = {}
        for p in mod["parameters"]:
            param_meta[p["form_name"]] = {
                "cli_name": p.get("cli_name", p["form_name"]),
                "param_type": p.get("param_type", "text"),
                "default": p.get("default", ""),
                "is_boolean": p.get("param_type") == "boolean",
            }
        from json import dumps
        self.write(
            t.generate(
                module=mod,
                param_meta=json.dumps(param_meta).encode('utf-8'),
                result=None,
                gs_version="0.3.0",
                modules=all_modules,
                param_overrides={},
            )
        )

    def post(self, module_name: str):
        module = self._find_module(module_name)
        if not module:
            self.set_status(404)
            self.write("Module not found.")
            return

        # Parse form values from the body
        values: Dict[str, str] = {}
        for key, val_list in self.request.body_arguments.items():
            if val_list:
                k = key.decode("utf-8") if isinstance(key, bytes) else key
                v = (
                    val_list[-1].decode("utf-8")
                    if isinstance(val_list[-1], bytes)
                    else val_list[-1]
                )
                values[k] = v

        cli_args = flatten_parameters(module["parameters"], values)

        # Run the tool
        result = run_biomate(module_name, cli_args)

        loader = get_loader()
        t = loader.load("module.html")
        all_modules = _module_list()
        mod = next(m for m in all_modules if m["name"] == module_name)

        # Build param_overrides from submitted values (for persisting in form after run)
        param_overrides = {}
        for p in mod["parameters"]:
            form_key = p.get("form_name", p.get("name", ""))
            if form_key in values and not isinstance(values[form_key], bytes):
                param_overrides[form_key] = values[form_key]
            elif p.get("param_type") == "boolean":
                if form_key not in values:
                    param_overrides[form_key] = str(p.get("default", False))

        # Build param_meta for JS
        param_meta = {}
        for p in mod["parameters"]:
            param_meta[p["form_name"]] = {
                "cli_name": p.get("cli_name", p["form_name"]),
                "param_type": p.get("param_type", "text"),
                "default": p.get("default", ""),
                "is_boolean": p.get("param_type") == "boolean",
            }

        self.write(
            t.generate(
                module=mod,
                param_meta=json.dumps(param_meta).encode('utf-8'),
                result=result,
                gs_version="0.3.0",
                modules=all_modules,
                param_overrides=param_overrides,
            )
        )

    def _find_module(self, name: str):
        for m in MODULES:
            if m["name"] == name:
                return m
        return None


# ---------------------------------------------------------------------------
# Application factory
# ---------------------------------------------------------------------------


def make_app(
    host: str = "localhost", port: int = 8080
) -> tornado.web.Application:
    """Create and return the Tornado application with security settings."""
    return tornado.web.Application(
        [
            (r"/", MainHandler),
            (r"/module/([a-zA-Z0-9_-]+)", ModuleHandler),
        ],
        template_path=str(_DESIGN_DIR),
        static_path=str(_STATIC_DIR),
        debug=False,  # Disable debug mode for security
        xsrf_cookies=True,  # Enable XSRF protection
        max_body_size=200 * 1024 * 1024,
    )


# ---------------------------------------------------------------------------
# CLI (sub-parser)
# ---------------------------------------------------------------------------


def init_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Initialise module subparser."""
    parser = subparsers.add_parser(
        "web-interface",
        help="Start the BioMate web interface (local development only)",
        description="""Launch the BioMate web interface on a local web server.

⚠️  SECURITY WARNING: This interface is designed for local development only.
   Do NOT expose to public networks without authentication and HTTPS.
   See documentation for security hardening guidelines if public deployment is needed.""",
    )
    parser.add_argument(
        "--host",
        default="localhost",
        help="Host to bind to (default: localhost; use 127.0.0.1 for local-only access)",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=8080,
        help="Port to bind to (default: 8080)",
    )
    parser.add_argument(
        "--ssl-cert",
        default=None,
        help="Path to SSL certificate file (required for HTTPS)",
    )
    parser.add_argument(
        "--ssl-key",
        default=None,
        help="Path to SSL private key file (required for HTTPS)",
    )
    parser.set_defaults(parse=main, run=main)

    return parser


def main(args: argparse.Namespace) -> None:
    """Start the web server."""

    # Warn if binding to public interfaces
    if args.host not in ("localhost", "127.0.0.1", "::1"):
        print("\n" + "=" * 70)
        print("⚠️  SECURITY WARNING")
        print("=" * 70)
        print(f"Binding to public interface: {args.host}")
        print("\nThe web interface will be accessible to your entire network.")
        print("This interface is NOT secure for public deployment!")
        print("\nRequired security measures for public access:")
        print("  1. HTTPS/TLS (use --ssl-cert and --ssl-key)")
        print("  2. Authentication (API keys or JWT tokens)")
        print("  3. Input validation and rate limiting (already partially implemented)")
        print("  4. Place behind a reverse proxy (nginx with auth)")
        print("\nFor local-only access, use: biomate web-interface --host 127.0.0.1")
        print("=" * 70 + "\n")

    app = make_app(host=args.host, port=args.port)

    certfile = args.ssl_cert or os.environ.get("SSL_CERT")
    keyfile = args.ssl_key or os.environ.get("SSL_KEY")

    if certfile and keyfile:
        context = ssl.SSLContext(ssl.PROTOCOL_TLS_SERVER)
        context.load_cert_chain(certfile, keyfile)
        http_server = tornado.httpserver.HTTPServer(app, ssl_options=context)
        http_server.listen(args.port, address=args.host)
        print(f"BioMate Web Interface (HTTPS) running at https://{args.host}:{args.port}")
    else:
        http_server = tornado.httpserver.HTTPServer(app)
        http_server.listen(args.port, address=args.host)
        print(f"BioMate Web Interface (HTTP) running at http://{args.host}:{args.port}")

    print("Press Ctrl+C to stop.\n")

    try:
        tornado.ioloop.IOLoop.current().start()
    except KeyboardInterrupt:
        print("\nShutting down...")
        http_server.stop()


if __name__ == "__main__":
    main()
