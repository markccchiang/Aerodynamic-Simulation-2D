#!/bin/bash
# Open the built Sphinx docs in the default browser, on macOS or Linux.
set -euo pipefail

# Resolve relative to this script, so it works from any directory.
index="$(cd "$(dirname "$0")" && pwd)/docs/_build/html/index.html"

if [ ! -f "$index" ]; then
    echo "Docs are not built yet. Build them first with:" >&2
    echo "  uv run sphinx-build -M html docs docs/_build" >&2
    exit 1
fi

case "$(uname -s)" in
    Darwin)
        open "$index"
        ;;
    Linux)
        if command -v xdg-open >/dev/null 2>&1; then
            xdg-open "$index"
        else
            echo "xdg-open not found; open this file in a browser:" >&2
            echo "  $index" >&2
            exit 1
        fi
        ;;
    *)
        echo "Unsupported OS; open this file in a browser:" >&2
        echo "  $index" >&2
        exit 1
        ;;
esac
