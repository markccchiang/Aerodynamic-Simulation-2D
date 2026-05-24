"""Entry point: launch the aerosim web UI.

    uv run python src/web.py                 # serve http://127.0.0.1:8000
    uv run python src/web.py --port 9000     # pick a port
"""

from __future__ import annotations

import argparse

import uvicorn

from aerosim.webapp import app


def main():
    ap = argparse.ArgumentParser(description="aerosim web UI server")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8000)
    args = ap.parse_args()
    print(f"aerosim web UI -> http://{args.host}:{args.port}")
    uvicorn.run(app, host=args.host, port=args.port)


if __name__ == "__main__":
    main()
