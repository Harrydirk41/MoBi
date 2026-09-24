"""PBPK Copilot: a local web app over the PBPK modeling agent.

Serves a browser UI and streams the real agent loop (inspect -> sweep ->
optimize -> grade) live over Server-Sent Events. Runs on the user's machine
because it needs PKSim.CLI and an Anthropic key; the browser talks only to
this local server.

    pip install -e ".[copilot]"
    set ANTHROPIC_API_KEY=...
    set PKPD_PKSIM_CLI=...\\PKSim.CLI.exe
    python -m copilot.server        # then open http://127.0.0.1:8765
"""
