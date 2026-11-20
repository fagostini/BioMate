# Web Interface

This sub-command starts a small local web server with one form per module. Each form builds the command line and runs `biomate <module> ...` on the server. The page shows the output and the exit code once the run finishes. Nspector HTML reports produced by a run can be opened from the result page.

> _**Warning:** The interface is meant for local use only. It has no authentication. Do not expose it to a network you do not trust. A host other than `localhost`, `127.0.0.1` or `::1` prints a security warning at start-up._

## Usage

The minimal command is `biomate web-interface`, which serves the interface at `http://localhost:8080`. Stop it with `Ctrl+C`.

## Options

- `--host`: Host to bind to. Default is `localhost`.
- `--port`: Port to bind to. Default is `8080`.
- `--ssl-cert`, `--ssl-key`: Certificate and private key files. When both are given, the server uses HTTPS. They take priority over the `SSL_CERT` and `SSL_KEY` environment variables, which are used when the options are not given.

## Notes

- The modules available in the interface are Blabber, Dirstruct, FastRewind, Index, Nspector, Strainer, Waif and Triage.
- Runs have no timeout: a long run keeps the request open until it finishes. For long runs, use the command line, where the progress is visible.
- Paths in the forms are paths on the machine that runs the server.
